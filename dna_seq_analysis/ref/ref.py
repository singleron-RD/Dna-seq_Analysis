import sys
from pathlib import Path
import subprocess
import tempfile
import unittest
from urllib.request import urlretrieve
from zipfile import ZipFile
from tempfile import NamedTemporaryFile

from dna_seq_analysis.tools.common import *



class Get_ref():
    def __init__(self,args, display_title=None):
        
        Step.__init__(self, args, display_title=display_title)
        
        self.config_path = args.config_path
        self.config = parse_config(self.config_path)
        self.species = args.species.lower() if args.species else self.config['ref']['species']
        self.build = args.build if args.build else self.config['ref']['build']
        self.release = args.release if args.release else self.config['ref']['release']
        self.outdir = Path(self.config['outdir'])
        
        # resources
        self.resources_dir = Path(self.config['genomedir'])
        self.resources_dir.mkdir(parents=True,exist_ok=True)
        # log
        self.log_dir = self.outdir/"logs/ref"
        self.log_dir.mkdir(parents=True,exist_ok=True)

        

    def get_genome(self):
        """
        curl genome
        """
        genome_fasta = self.resources_dir/"genome.fasta"
        genome_fasta_fai = self.resources_dir/"genome.fasta.fai"
        genome_fasta_dict = self.resources_dir/"genome.dict"
        genome_gtf = self.resources_dir/"genome.gtf"
        genome_refFlat = self.resources_dir/"genome.refFlat"


        branch = ""
        if self.release >= 81 and self.build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"

        spec = (f"{self.build}" if int(self.release) > 75 else f"{self.build}.{self.release}")

        suffixes = ""
        datatype = self.args.datatype
        chromosome = self.args.chromosome
        
        if datatype == "dna":
            if chromosome:
                suffixes = [f"dna.chromosome.{chromosome}.fa.gz"]
            else:
                suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
        elif datatype == "cdna":
            suffixes = ["cdna.all.fa.gz"]
        elif datatype == "cds":
            suffixes = ["cds.all.fa.gz"]
        elif datatype == "ncrna":
            suffixes = ["ncrna.fa.gz"]
        elif datatype == "pep":
            suffixes = ["pep.all.fa.gz"]
        else:
            raise ValueError("this datatype is invalid, must be one of dna, cdna, cds, ncrna, pep")

        if chromosome:
            if not datatype == "dna":
                raise ValueError(
                    "invalid datatype, to select a single chromosome the datatype must be dna"
                )

        success = False
        for suffix in suffixes:
            url = f"ftp://ftp.ensembl.org/pub/{branch}release-{self.release}/fasta/{self.species}/{datatype}/{self.species.capitalize()}.{spec}.{suffix}"
            try:
                _cmd = (f"curl -sSf {url} > /dev/null 2> /dev/null")
                debug_subprocess_call(_cmd)     
            except subprocess.CalledProcessError:
                continue
            cmd = (f"(curl -L {url} | gzip -d > {genome_fasta})")
            if not Path({genome_fasta}).exists():
                debug_subprocess_call(cmd)
            success = True
            break
        
        url = f"ftp://ftp.ensembl.org/pub/{branch}release-{self.release}/gtf/{self.species}/{self.species.capitalize()}.{self.build}.{self.release}.gtf.gz"
        cmd = (f"(curl -L {url} |gzip -d > {genome_gtf})")
        if not Path({genome_gtf}).exists():
            debug_subprocess_call(cmd)         
        cmd_gtf2refflat = (f"gtfToGenePred -genePredExt -geneNameAsName2 {genome_gtf} -ignoreGroupsWithoutExons /dev/stdout| "
                           r'''awk 'BEGIN { OFS="\t" } {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > '''
                           f"{genome_refFlat}")
        debug_subprocess_call(cmd_gtf2refflat) 
        

        if not success:
            print(
                "Unable to download requested sequence data from Ensembl. "
                "Did you check that this combination of species, build, and release is actually provided?",
                file=sys.stderr,
            )
            exit(1)

        faidx_dict = (f"samtools faidx {genome_fasta} > {genome_fasta_fai};samtools dict {genome_fasta} > {genome_fasta_dict}")
        debug_subprocess_call(faidx_dict)


    def get_known_variation(self):
        """
        get_known_variation
        """
        # input
        genome_fasta_fai = self.resources_dir/"genome.fasta.fai"
        # output
        variation = self.resources_dir/"variation.vcf.gz"
        log = self.log_dir/"get-known-variants.log"

        _type = self.args.type
        if self.release < 98:
            raise ValueError("Ensembl releases <98 are unsupported.")
    
        branch = ""
        if self.release >= 81 and self.build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"

        if _type == "all":
            if self.species == "homo_sapiens" and self.release >= 93:
                suffixes = [f"-chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y", "MT"]]
            else:
                suffixes = [""]
        elif _type == "somatic":
            suffixes = ["_somatic"]
        elif _type == "structural_variations":
            suffixes = ["_structural_variations"]
        else:
            raise ValueError(f"Unsupported type {_type} (only all, somatic, structural_variations are allowed)")

        species_filename = self.species if self.release >= 91 else self.species.capitalize()

        urls = [
            f"ftp://ftp.ensembl.org/pub/{branch}release-{self.release}/variation/vcf/{self.species}/{species_filename}{suffix}.{ext}"
            for suffix in suffixes
            for ext in ["vcf.gz", "vcf.gz.csi"]
        ]
        names = [Path(url).name for url in urls if url.endswith(".gz")]

        try:  
            gather = "curl {urls}".format(urls=" ".join(map("-O {}".format, urls)))
            with tempfile.TemporaryDirectory() as tmpdir:
                if genome_fasta_fai.exists():
                    cmd = (
                        f"(cd {tmpdir};{gather} && "  
                        f"bcftools concat -Oz --naive {' '.join(names)} > concat.vcf.gz && "
                        f"bcftools reheader --fai {genome_fasta_fai} concat.vcf.gz "
                        f"> {variation})"
                    )
                else:
                    cmd = (
                        f"(cd {tmpdir};{gather} && "
                        f"bcftools concat -Oz --naive {' '.join(names)} "
                        f"> {variation})"
                    )
                debug_subprocess_call(cmd)
        except subprocess.CalledProcessError as e:
            if log:
                sys.stderr = open(log, "a")
            print(
                "Unable to download variation data from Ensembl. "
                "Did you check that this combination of species, build, and release is actually provided? ",
                file=sys.stderr,
            )
            exit(1)

    def remove_iupac_codes(self):
        # input
        variation = self.resources_dir/"variation.vcf.gz"
        # output
        noiupac_variation = self.resources_dir/"variation.noiupac.vcf.gz"
        log = self.log_dir/"tabix-variation.log"

        cmd = (f"rbt vcf-fix-iupac-alleles < {variation} |bcftools view -Oz > {noiupac_variation};tabix -p vcf {noiupac_variation} 2>{log}")
        debug_subprocess_call(cmd)

    def bwa_index(self):
        genome_fasta = self.resources_dir/"genome.fasta"
        
        # Prefix that should be used for the database
        prefix = self.args.index_prefix if self.args.index_prefix else ''

        if prefix:
            prefix = "-p " + prefix

        # Contrunction algorithm that will be used to build the database, default is bwtsw
        construction_algorithm = self.args.algorithm if self.args.algorithm else ''

        if construction_algorithm:
            construction_algorithm = "-a " + construction_algorithm

        cmd = ( "bwa index " 
                f"{prefix} " 
                f"{construction_algorithm} " 
                f"{genome_fasta} " )
        debug_subprocess_call(cmd)


    def get_vep_cache(self):
        # output
        vep_cache_dir = self.resources_dir/"vep/cache"
        vep_cache_dir.mkdir(parents=True,exist_ok=True)
        
        branch = ""
        if self.release >= 81 and self.build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"
        species_filename = self.species if self.release >= 91 else self.species.capitalize()
        vep_cache =  'indexed-vep_cache' if self.release >= 91 else 'VEP'
         
        url = f"ftp://ftp.ensembl.org/pub/{branch}release-{self.release}/variation/{vep_cache}/{species_filename}_vep_{self.release}_{self.build}.tar.gz"
        cmd = (f'curl -L {url} -o {vep_cache_dir}/{species_filename}_vep_{self.release}_{self.build}.tar.gz &&'
               f'cd {vep_cache_dir} && tar -xzf {species_filename}_vep_{self.release}_{self.build}.tar.gz && rm {species_filename}_vep_{self.release}_{self.build}.tar.gz')
        if not Path(f'{vep_cache_dir}/{species_filename}_vep_{self.release}_{self.build}.tar.gz').exists():
            debug_subprocess_call(cmd)


    def get_vep_plugins(self):
        print('vep_plugins')
        # output
        vep_plugins_dir = self.resources_dir/"vep/plugins"
        vep_plugins_dir.mkdir(parents=True,exist_ok=True)
        # log
        log = self.log_dir/"vep-plugins.log"

        if log:
            sys.stderr = open(log, "w")

        with NamedTemporaryFile() as tmp:
            urlretrieve(f"https://github.com/Ensembl/VEP_plugins/archive/release/{self.release}.zip",tmp.name,)
            with ZipFile(tmp.name) as f:
                for member in f.infolist():
                    memberpath = Path(member.filename)
                    if len(memberpath.parts) == 1:
                        # skip root dir
                        continue
                    targetpath = vep_plugins_dir / memberpath.relative_to(memberpath.parts[0])
                    if member.is_dir():
                        targetpath.mkdir()
                    else:
                        with open(targetpath, "wb") as out:
                            out.write(f.read(member.filename))


    @add_log
    def main(self):
        self.get_genome()
        self.get_known_variation()
        self.remove_iupac_codes()
        self.bwa_index()
        self.get_vep_cache()
        self.get_vep_plugins()


@add_log
def ref(args):
    runner = Get_ref(args)
    runner.main()
    
    
def get_opts_ref(parser, sub_program=True):
    parser.add_argument('--species',help="Ensembl species name.",required=True)
    parser.add_argument('--release',help="Ensembl release.",required=True,type=int)
    parser.add_argument('--build',help="Genome build.",required=True)
    parser.add_argument('--datatype',help="Sequence types.",choices=['dna','cdna','cds','ncrna','pep'],default='dna')
    parser.add_argument('--chromosome',help="Select a specific chromosome for analysis.")
    parser.add_argument('--type',help="Ensembl VCF (Variant Call Format) files types.",choices=['all','somatic','structural_variations'],default='all')
    parser.add_argument('-p','--index_prefix',help="The prefix of the generated index files.(same as fasta name)")
    parser.add_argument('-a','--algorithm',help="BWT construction algorithm: bwtsw, is or rb2 [auto]")
    
    if sub_program:
        parser = s_common(parser)
    return parser


if __name__ == '__main__':
    unittest.main()
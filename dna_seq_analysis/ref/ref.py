import sys
from pathlib import Path
import subprocess
import os
import tempfile
from urllib.request import urlretrieve
from zipfile import ZipFile
from tempfile import NamedTemporaryFile

from dna_seq_analysis.tools.common import *



class Get_ref():
    def __init__(self):
        self.species = config['ref']['species'].lower()
        self.build = config['ref']['build']
        self.release = config['ref']['release']
        self.outdir = config['outdir']

        # resources
        self.resources_dir = resource_dir
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


        branch = ""
        if self.release >= 81 and self.build == "GRCh37":
            # use the special grch37 branch for new releases
            branch = "grch37/"

        spec = (f"{self.build}" if int(self.release) > 75 else f"{self.build}.{self.release}")

        suffixes = ""
        datatype = config['ref'].get("datatype", "") if config['ref'].get("datatype", "") else "dna"
        chromosome = config['ref'].get("chromosome", "")
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
            debug_subprocess_call(cmd)
            success = True
            break

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
        genome_fai = self.resources_dir/"genome.fasta.fai"
        # output
        variation = self.resources_dir/"variation.vcf.gz"
        log = self.log_dir/"get-known-variants.log"

        _type = config['ref'].get("type", "") if config['ref'].get("type", "") else "all"
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
        names = [os.path.basename(url) for url in urls if url.endswith(".gz")]

        try:  
            gather = "curl {urls}".format(urls=" ".join(map("-O {}".format, urls)))
            workdir = os.getcwd()
            with tempfile.TemporaryDirectory() as tmpdir:
                if genome_fai:
                    cmd = (
                        f"(cd {tmpdir};{gather} && "  
                        f"bcftools concat -Oz --naive {names} > concat.vcf.gz && "
                        f"bcftools reheader --fai {workdir}/{genome_fai} concat.vcf.gz "
                        f"> {workdir}/{variation})"
                    )
                else:
                    cmd = (
                        f"(cd {tmpdir};{gather} && "
                        f"bcftools concat -Oz --naive {names} "
                        f"> {workdir}/{variation})"
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
        
        # Check inputs/arguments.
        if len(genome_fasta) == 0:
            raise ValueError("A reference genome has to be provided!")
        elif len(genome_fasta) > 1:
            raise ValueError("Only one reference genome can be inputed!")

        # Prefix that should be used for the database
        prefix = config['ref'].get("prefix", "")

        if len(prefix) > 0:
            prefix = "-p " + prefix

        # Contrunction algorithm that will be used to build the database, default is bwtsw
        construction_algorithm = config['ref'].get("algorithm", "")

        if len(construction_algorithm) != 0:
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

        cmd = (
            "vep_install --AUTO cf "
            f"--SPECIES {self.species} "
            f"--ASSEMBLY {self.build} "
            f"--VERSION {self.release} "
            f"--CACHEDIR {str(vep_cache_dir)} "
            "--CONVERT "
            "--NO_UPDATE "
            )
        debug_subprocess_call(cmd)


    def get_vep_plugins(self):
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


runner = Get_ref()
runner.main()
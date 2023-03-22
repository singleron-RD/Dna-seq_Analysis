from collections import namedtuple
import numbers
import subprocess
import argparse



class Picard:
    def __init__(self,args):
        self.bam = args.bam
        self.refflat = args.refflat
        self.sample = args.sample
        self.outdir = args.outdir

        self.Metric = namedtuple("Metric", "name value total fraction")
        self.metric_list = []
        self.picard_region_log = f'{self.outdir}/{self.sample}_region.log'
        self.stat_file = f'{self.outdir}/{self.sample}_stat.txt'

    def picard(self):
        cmd = [
            'picard',
            'CollectRnaSeqMetrics',
            '-XX:ParallelGCThreads=6',
            '-Xmx20G',
            f'I={self.bam}',
            f'O={self.picard_region_log}',
            f'REF_FLAT={self.refflat}',
            'STRAND=FIRST_READ_TRANSCRIPTION_STRAND'
        ]
        cmd_str = ' '.join(cmd)
        subprocess.check_call(cmd_str,shell=True)
        
    def add_metric(self,name, value=None, total=None, fraction=None):
            '''add metric to metric_list
            '''
            self.metric_list.append(self.Metric(
                name=name, value=value, total=total, fraction=fraction
            ))

    def add_other_metrics(self):
            """
            add picard region bases
            add region plot
            if debug, add ribosomal RNA reads percent
            """

            with open(self.picard_region_log, 'r') as picard_log:
                region_dict = {}
                for line in picard_log:
                    if not line:
                        break
                    if line.startswith('## METRICS CLASS'):
                        header = picard_log.readline().strip().split('\t')
                        data = picard_log.readline().strip().split('\t')
                        region_dict = dict(zip(header, data))
                        break

            total = float(region_dict['PF_ALIGNED_BASES'])
            exonic_regions = int(region_dict['UTR_BASES']) + int(region_dict['CODING_BASES'])
            exonic_regions_fraction = float(exonic_regions/total)
            intronic_regions = int(region_dict['INTRONIC_BASES'])
            intronic_regions_fraction = float(intronic_regions/total)
            intergenic_regions = int(region_dict['INTERGENIC_BASES'])
            intergenic_regions_fraction = float(intergenic_regions/total)

            self.add_metric(
                name='Base Pairs Mapped to Exonic Regions',
                value=exonic_regions,
                total=total,
                fraction=exonic_regions_fraction
            )
            self.add_metric(
                name='Base Pairs Mapped to Intronic Regions',
                value=intronic_regions,
                total=total,
                fraction=intronic_regions_fraction
            )
            self.add_metric(
                name='Base Pairs Mapped to Intergenic Regions',
                value=intergenic_regions,
                total=total,
                fraction=intergenic_regions_fraction
            )

    def metric_list_to_stat(self):
        with open(self.stat_file, 'w') as stat_handle:
            for metric in self.metric_list:
                line = f'{metric.name}: '
                value = metric.value
                fraction = metric.fraction
                value_bool = value == 0 or value
                fraction_bool = fraction == 0 or fraction
                if fraction_bool:
                    fraction = round(fraction * 100, 2)
                if value_bool:
                    if isinstance(value, numbers.Number):
                        line += format(value, ',')
                        if fraction_bool:
                            line += f'({fraction}%)'
                    else:
                        line += value
                elif fraction_bool:
                    line += f'{fraction}%'
                stat_handle.write(line + '\n')

    def run(self):
        self.picard()
        self.add_other_metrics()
        self.metric_list_to_stat()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='picard_args')
    parser.add_argument('--bam',help='mapped bam', required=True)
    parser.add_argument('--outdir', help='Output diretory.', required=True)
    parser.add_argument('--refflat',help='ref genome', required=True)
    parser.add_argument('--sample',help='sample name', required=True)
    args = parser.parse_args()
    run_picard = Picard(args)
    run_picard.run()
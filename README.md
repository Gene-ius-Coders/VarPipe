# VarPipe

VarPipe is a pipeline that utilizes VarScan to detect variants and discover somatic mutations/CNVs in next-generation sequencing data. **(alignment, varscan, and trimming)**

## Installation

Put all of the files in a folder. Then, from [here](https://github.com/arq5x/bedtools2/releases), download the bedtools binary file and place it in the folder.
## Usage
```bash
./pipeline.sh -i [data_location.fastq] -r [reference_location.fasta] -o [output_location] -t true -b true -p 0.05
```
More examples can be found as follows:
```bash
./pipeline.sh -h
```


## Contributing

Pull requests are welcome.

## License
Project is distributed under [MIT](https://choosealicense.com/licenses/mit/) licence.

#! /bin/bash
#arguments
while test $# -gt 0; do
           case "$1" in
                -h|--help)
                    echo "T group pipeline **(alignment, varscan, trimming)**"
                    echo " "
                    echo "options:"
                    echo "      -h, --help                    show brief help"
                    echo "      -r, --reference               specify reference"
                    echo "      -i, --input                   specify input folder"
                    echo "      -o, --output-dir              specify a directory to store output in"
                    echo "      -p,                           specify p_value (optional - default=0.01)"
                    echo "      -t,                           trimming (optional true or false - default=false)"
                    echo "      -b,                           bed file (optional true or false - default=false)"
                    echo "examples:"
                    echo "      1| ./pipeline.sh -i ./data -r reference/chr21_22_v37.fasta -o ./output"
                    echo "      2| ./pipeline.sh -i ./data -r reference/chr21_22_v37.fasta -o ./output -t true"
                    echo "      3| ./pipeline.sh -i ./data -r reference/chr21_22_v37.fasta -o ./output -t true -p 0.05"
                    echo "      4| ./pipeline.sh -i ./data -r reference/chr21_22_v37.fasta -o ./output -t true -b true -p 0.05"
                    echo "note:"
                    echo "      **skewer and bedtools binaries must be in the same path as pipeline.sh, same for picard.js and varscan.js**"
                    exit 0
                    ;;
                -i|--input)
                    shift
                    if test $# -gt 0; then
                        cases_folder=$1
                    else
                        echo "no input folder specified"
                        exit 1
                    fi
                    shift
                    ;;
                -r|--reference)
                    shift
                    if test $# -gt 0; then
                        reference_file=$1
                    else
                        echo "no reference file specified"
                        exit 1
                    fi
                    shift
                    ;;
                -o|--output_dir)
                    shift
                    if test $# -gt 0; then
                       output_dir=$1
                    else
                        echo "no output directory specified"
                        exit 1
                    fi
                    shift
                    ;;
                -p)
                    shift
                    p_value=$1
                    shift
                    ;;
                -t)
                    shift
                    trimming=$1
                    shift
                    ;;
                -b)
                    shift
                    bedtools=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   exit 1
                   ;;
          esac
  done
if [ -z "$cases_folder" ]
then
    echo "input folder is required"
    exit 1
fi
if [ -z "$reference_file" ]
then
    echo "reference file is required"
    exit 1
fi
if [ -z "$output_dir" ]
then
    echo "output directory is required"
    exit 1
fi

#set default p-val
if [ -z "$p_value" ]
then
    p_value=0.01
fi

#trimming permission
if [ -z "$trimming" ]
then
    trimming="false"
else
    SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
    sudo chmod a+x "$SCRIPTPATH"/skewer-0.2.2-linux-x86_64
    skewer_path="$SCRIPTPATH"/skewer-0.2.2-linux-x86_64
fi

#bedtools permission
if [ -z "$bedtools" ]
then
    bedtools="false"
else
    SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
    sudo chmod a+x "$SCRIPTPATH"/bedtools
    bedtools="$SCRIPTPATH"/bedtools
fi

#varscan binary path
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
sudo chmod a+x "$SCRIPTPATH"/VarScan.v2.3.9.jar
varscan="$SCRIPTPATH"/VarScan.v2.3.9.jar

#folder name of cases
name_of_cases=`for dir in "$cases_folder"/*; do basename "$dir"; done`
set -o noglob
IFS=$'\n' name_of_cases_arr=($name_of_cases)
set +o noglob

##color for echo
RED='\033[0;31m'
BGreen='\033[1;32m'
NC='\033[0m' # No Color

for name_of_a_case in "${name_of_cases_arr[@]}"  
do
    #make output directory
    mkdir -p "$output_dir/$name_of_a_case"

    #names of two fastq files
    paired_cases=`for dir in "$cases_folder/$name_of_a_case"/*.fastq*; do basename "$dir"; done`
    set -o noglob
    IFS=$'\n' paired_cases_arr=($paired_cases)
    set +o noglob

    #run fastqc for first file
    fastqc "$cases_folder/$name_of_a_case/${paired_cases_arr[0]}" --extract -o "$output_dir/$name_of_a_case/"
    INPUT1="${paired_cases_arr[0]}"
    SUBSTRING1=$(echo $INPUT1| cut -d'.' -f 1)
    flag1=`cat "$output_dir/$name_of_a_case/$SUBSTRING1""_fastqc/summary.txt" | awk '{printf ($2=="Adapter" && $1=="FAIL") ? $1"\n" : "";}'`
    err1=`cat "$output_dir/$name_of_a_case/$SUBSTRING1""_fastqc/summary.txt" | awk '{printf ($1=="FAIL") ? $0"\n" : "";}'`
    if [ "$err1" != "" ]
    then
        echo -e "${RED}((!WARNING)) for data: $cases_folder/$name_of_a_case/${paired_cases_arr[0]} ${NC}"
        echo -e "${RED} $err1 ${NC}"
        read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
    fi

    #run fastqc for second file
    fastqc "$cases_folder/$name_of_a_case/${paired_cases_arr[1]}" --extract -o "$output_dir/$name_of_a_case/"
    INPUT2="${paired_cases_arr[1]}"
    SUBSTRING2=$(echo $INPUT2| cut -d'.' -f 1)
    flag2=`cat "$output_dir/$name_of_a_case/$SUBSTRING2""_fastqc/summary.txt" | awk '{printf ($2=="Adapter" && $1=="FAIL") ? $1"\n" : "";}'`
    err2=`cat "$output_dir/$name_of_a_case/$SUBSTRING2""_fastqc/summary.txt" | awk '{printf ($1=="FAIL") ? $0"\n" : "";}'`
    if [ "$err2" != "" ]
    then
        echo -e "${RED}((!WARNING)) for data: $cases_folder/$name_of_a_case/${paired_cases_arr[1]} ${NC}"
        echo -e "${RED} $err2 ${NC}"
        read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
    fi

    #trimming
    if [ "$trimming" != "false" ] && [ "$flag1" == "FAIL" ] && [ "$flag2" == "FAIL" ]
    then
        trimming=`awk -v lengths=0 -v seq=0 'NR%4 == 2 {lengths+=length($0)} NR%4 == 2 {seq++} END {print lengths / seq}' "$cases_folder/$name_of_a_case/${paired_cases_arr[0]}"`
        $skewer_path -x AGATCGGAAGAG -y AGATCGGAAGAG -m pe -l "$trimming" "$cases_folder/$name_of_a_case/${paired_cases_arr[0]}" "$cases_folder/$name_of_a_case/${paired_cases_arr[1]}" 
    fi
done

for name_of_a_case in "${name_of_cases_arr[@]}"  
do
    #read trimmed data else main data
    paired_cases=`for dir in "$cases_folder/$name_of_a_case"/*-trimmed-pair*.fastq*; do basename "$dir"; done`
    if [ "$paired_cases" == "*-trimmed-pair*.fastq*" ]
    then
        paired_cases=`for dir in "$cases_folder/$name_of_a_case"/*.fastq*; do basename "$dir"; done`
    fi
    set -o noglob
    IFS=$'\n' paired_cases_arr=($paired_cases)
    set +o noglob
    bwa mem "$reference_file" "$cases_folder/$name_of_a_case/${paired_cases_arr[0]}" "$cases_folder/$name_of_a_case/${paired_cases_arr[1]}" > "$output_dir/$name_of_a_case/$name_of_a_case".untrimmed.sam
    samtools view -bT "$reference_file" "$output_dir/$name_of_a_case/$name_of_a_case".untrimmed.sam -o "$output_dir/$name_of_a_case/$name_of_a_case".bam
    
    #quality check of bam file
    validation=`java -jar ./picard.jar ValidateSamFile -I "$output_dir/$name_of_a_case/$name_of_a_case".bam -MODE SUMMARY --QUIET true | awk '{printf(/ERROR/) ? $0 : ""}'`
    if [ "$validation" != "" ]
    then
        echo -e "${RED} $validation ${NC}"
        read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
    fi

    #sort, index
    samtools sort "$output_dir/$name_of_a_case/$name_of_a_case".bam -o "$output_dir/$name_of_a_case/$name_of_a_case".sorted.bam
    samtools index "$output_dir/$name_of_a_case/$name_of_a_case".sorted.bam
    samtools idxstats "$output_dir/$name_of_a_case/$name_of_a_case".sorted.bam
    
    #varscan
    samtools mpileup -f "$reference_file" "$output_dir/$name_of_a_case/$name_of_a_case".sorted.bam | java -jar "$varscan" mpileup2cns --p-value "$p_value" --output-vcf 1 --variants 0 > "$output_dir/$name_of_a_case/$name_of_a_case".vcf

    #bed
    bedfile=`for dir in "$cases_folder/$name_of_a_case"/*.bed*; do basename "$dir"; done`
    if [ "$bedfile" != "*.bed*" ]
    then
        ./bedtools intersect -a "$output_dir/$name_of_a_case/$name_of_a_case".vcf -b "$cases_folder/$name_of_a_case/${bedfile[0]}" -header > "$output_dir/$name_of_a_case/$name_of_a_case"_filtered.vcf
    fi
done
echo "${BGreen}For annotation please go to 'https://wannovar.wglab.org/' and use the _filtered.vcf file, May God bless you! ${NC}"
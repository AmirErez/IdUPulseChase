ls -l *.txt | perl -ne 's/.*\s(\d+)hr_r(\d+)_m(\d+)_(\d+….*)/mv $1hr_r$2_m$3_$4 Ly6G+_$1hr_r$2_m$3.csv/og;print;'


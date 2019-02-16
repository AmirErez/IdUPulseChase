ls -l PreToImmatureB*.csv | perl -ne 's/.*PreToImmatureB_(\d+)hr_r(\d+)_m(\d+)\.csv/PreToImmatureB_$1hr_r$2_m$3 $1 $2 BM_PreToImmatureB/og;print;'


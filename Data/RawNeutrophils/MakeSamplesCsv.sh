ls -l Ly6G*.csv | perl -ne 's/.*Ly6G\+_(\d+)hr_r(\d+)_m(\d+)\.csv/Ly6G+_$1hr_r$2_m$3 $1 $2 BM_Ly6G+/og;print;'


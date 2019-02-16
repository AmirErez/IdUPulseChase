ls -l ImmTransMarg*.csv | perl -ne 's/.*ImmTransMarg_(\d+)hr_r(\d+)_m(\d+)\.csv/ImmTransMarg_$1hr_r$2_m$3 $1 $2 BM_ImmTransMarg/og;print;'


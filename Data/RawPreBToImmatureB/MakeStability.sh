cat rawStabilityPreToImmature.csv | perl -ne 's/\r/\n/og; s/\t/,/og; s/\d+:\s(.*m\d+).*,(\S*),(\S*),(\S*),(\S*)/PreToImmatureB_$1, $2, $3, $4, $5/og; print;' > StabilityPreToImmature.csv

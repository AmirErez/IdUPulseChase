cat rawStabilityNeut.csv | perl -ne 's/\r/\n/og; s/\t/,/og; s/\d+:\s(.*m\d+).*,(\S*),(\S*),(\S*),(\S*)/Ly6G+_$1, $2, $3, $4, $5/og; print;' > StabilityNeut.csv

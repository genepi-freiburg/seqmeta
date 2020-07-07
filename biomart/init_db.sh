tail -n+2 mart_export_exons.txt > mart_export_exons_noheader.txt
rm ensembl_exons.sqlite
sqlite3 ensembl_exons.sqlite < ensembl_exons.sql
sqlite3 ensembl_exons.sqlite ".separator \t" ".import mart_export_exons_noheader.txt exons"
rm mart_export_exons_noheader.txt

echo "MSU7.Chr4.mPing 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing > MSU7.Chr4.mPing.500.list
perl list2csv.pl --list MSU7.Chr4.mPing.500.list --project MSU7.Chr4.mPing.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --lib in_libs.MSU7.Chr4.mPing.500.csv --group in_groups.MSU7.Chr4.mPing.500.csv --project MSU7.Chr4.mPing.500
bash MSU7.Chr4.mPing.500.sh &
perl RedistributeBam.py --input .

echo "MSU7.Chr4.mPing 200bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing > MSU7.Chr4.mPing.200.list
perl list2csv.pl --list MSU7.Chr4.mPing.200.list --project MSU7.Chr4.mPing.200
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --lib in_libs.MSU7.Chr4.mPing.200.csv --group in_groups.MSU7.Chr4.mPing.200.csv --project MSU7.Chr4.mPing.200
bash MSU7.Chr4.mPing.200.sh &
perl RedistributeBam.py --input .

echo "MSU7.Chr4.ALL 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL > MSU7.Chr4.ALL.500.list
perl list2csv.pl --list MSU7.Chr4.ALL.500.list --project MSU7.Chr4.ALL.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --lib in_libs.MSU7.Chr4.ALL.500.csv --group in_groups.MSU7.Chr4.ALL.500.csv --project MSU7.Chr4.ALL.500
bash MSU7.Chr4.ALL.500.sh &
perl RedistributeBam.py --input .

echo "MSU7.Chr3.mPing 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.mPing > MSU7.Chr3.mPing.500.list
perl list2csv.pl --list MSU7.Chr3.mPing.500.list --project MSU7.Chr3.mPing.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --lib in_libs.MSU7.Chr3.mPing.500.csv --group in_groups.MSU7.Chr3.mPing.500.csv --project MSU7.Chr3.mPing.500
bash MSU7.Chr3.mPing.500.sh &
perl RedistributeBam.py --input .

echo "MSU7.Chr3.ALL 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL > MSU7.Chr3.ALL.500.list
perl list2csv.pl --list MSU7.Chr3.ALL.500.list --project MSU7.Chr3.ALL.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --lib in_libs.MSU7.Chr3.ALL.500.csv --group in_groups.MSU7.Chr3.ALL.500.csv --project MSU7.Chr3.ALL.500
bash MSU7.Chr3.ALL.500.sh &
perl RedistributeBam.py --input .


echo "FLY603.Chr2L.pogo 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/FLY603.Chr2L.pogo > FLY603.Chr2L.pogo.500.list
perl list2csv.pl --list FLY603.Chr2L.pogo.500.list --project FLY603.Chr2L.pogo.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/FLY603.Chr2L.fa --lib in_libs.FLY603.Chr2L.pogo.500.csv --group in_groups.FLY603.Chr2L.pogo.500.csv --project FLY603.Chr2L.pogo.500
bash FLY603.Chr2L.pogo.500.sh &
perl RedistributeBam.py --input .

echo "TAIR10.Chr1.mPing 500bp library"
python PrepareCSV.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/TAIR10.Chr1.mPing > TAIR10.Chr1.mPing.500.list
perl list2csv.pl --list TAIR10.Chr1.mPing.500.list --project TAIR10.Chr1.mPing.500
perl runMapping.pl --ref /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/TAIR10.Chr1.fa --lib in_libs.TAIR10.Chr1.mPing.500.csv --group in_groups.TAIR10.Chr1.mPing.500.csv --project TAIR10.Chr1.mPing.500
bash TAIR10.Chr1.mPing.500.sh &
perl RedistributeBam.py --input .

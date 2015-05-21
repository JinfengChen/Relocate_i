echo "MSU7_Chr4"
perl /rhome/cjinfeng/software/bin/fastaDeal.pl --pat Chr4 ./Reference/MSU_r7.fa > ./Reference/MSU7.Chr4.fa
perl /rhome/cjinfeng/software/bin/fastaDeal.pl --reform line50 MSU7.Chr4.fa > ./Reference/MSU7.Chr4a.fa
mv MSU7.Chr4a.fa MSU7.Chr4.fa
grep "Chr4" /rhome/cjinfeng/BigData/00.RD/GenomeAlign/Lastz/input/mask/MSU_r7.fa.RepeatMasker.out.bed > MSU7.Chr4.fa.RepeatMasker.out.bed
grep "Chr4" /rhome/cjinfeng/BigData/00.RD/GenomeAlign/Lastz/input/mask/MSU_r7.fa.RepeatMasker.out.existingTE > MSU7.Chr4.fa.RepeatMasker.out.existingTE

echo "MSU7_Chr3"
perl /rhome/cjinfeng/software/bin/fastaDeal.pl --pat Chr3 MSU_r7.fa > MSU7.Chr3.fa
perl /rhome/cjinfeng/software/bin/fastaDeal.pl --reform line50 MSU7.Chr3.fa > ./Reference/MSU7.Chr3a.fa
mv MSU7.Chr3a.fa MSU7.Chr3.fa
grep "Chr3" /rhome/cjinfeng/BigData/00.RD/GenomeAlign/Lastz/input/mask/MSU_r7.fa.RepeatMasker.out.bed > MSU7.Chr3.fa.RepeatMasker.out.bed

echo "Flybase dmel6.03 Chr2L"
perl /rhome/cjinfeng/software/bin/fastaDeal.pl --pat 2L dmel-all-chromosome-r6.03.fasta > FLY603.Chr2L.fa
cp /rhome/cjinfeng/BigData/00.RD/seqlib/Fly/FLY603.Chr2L.fa.RepeatMasker.out.bed ./
cp /rhome/cjinfeng/BigData/00.RD/seqlib/Fly/FLY603.Chr2L.fa.RepeatMasker.out ./

echo "TAIR10 Chr1"
cp /rhome/cjinfeng/BigData/00.RD/seqlib/Arabidopsis/TAIR10_chr1.fas TAIR10.Chr1.fa
grep "Chr1" /rhome/cjinfeng/BigData/00.RD/seqlib/Arabidopsis/TAIR10.fa.RepeatMasker.out.bed > TAIR10.Chr1.fa.RepeatMasker.out.bed


echo "MSU7_Chr4_Flybase_Chr2L mPing, pogo"
cat MSU7.Chr4.fa FLY603.Chr2L.fa > Mix.ALL.fa
cat mping.fa pogo.fa > repeat.fa
mkdir Mix.ALL.Repeat
cat ../MSU7.Chr4.mPing/MSU7.Chr4.mPing.rep1.fasta ../FLY603.Chr2L.pogo/FLY603.Chr2L.pogo.rep1.fasta > Mix.ALL.Repeat.rep1.fasta
cat ../MSU7.Chr4.mPing/MSU7.Chr4.mPing.rep1.gff ../FLY603.Chr2L.pogo/FLY603.Chr2L.pogo.rep1.gff > Mix.ALL.Repeat.rep1.gff

python Simulation_illumina.py --input ./Mix.ALL.Repeat --size 500


echo "MSU7_Chr4_mPing"
echo "Simulate genome"
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr4 --number 200 --prefix MSU7.Chr4.mPing.rep1
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr4 --number 200 --prefix MSU7.Chr4.mPing.rep2
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr4 --number 200 --prefix MSU7.Chr4.mPing.rep3

echo "Simulate illumina reads"
python Simulation_illumina.py --input ./MSU7.Chr4.mPing --size 500
python Simulation_illumina.py --input ./MSU7.Chr4.mPing --size 350
python Simulation_illumina.py --input ./MSU7.Chr4.mPing --size 200


echo "MSU7_Chr4_ALL"
echo "Simulate genome"
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr4 --number 200 --prefix MSU7.Chr4.ALL.rep1
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr4 --number 200 --prefix MSU7.Chr4.ALL.rep2
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr4 --number 200 --prefix MSU7.Chr4.ALL.rep3

echo "Simulate illumina reads"
python Simulation_illumina.py --input ./MSU7.Chr4.ALL --size 500


echo "MSU7_Chr3_mPing"
echo "Simulate genome"
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr3 --number 200 --prefix MSU7.Chr3.mPing.rep1
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr3 --number 200 --prefix MSU7.Chr3.mPing.rep2
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/MSU_r7.fa --te mPing  --chr Chr3 --number 200 --prefix MSU7.Chr3.mPing.rep3
 
echo "Simulate illumina reads"
python Simulation_illumina.py --input ./MSU7.Chr3.mPing --size 500


echo "MSU7_Chr3_ALL"
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr3 --number 200 --prefix MSU7.Chr3.ALL.rep1
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr3 --number 200 --prefix MSU7.Chr3.ALL.rep2
python Simulation_MultiTE_Insertion.py --repeat ./Reference/RiceTE.TSD.fa --genome ./Reference/MSU_r7.fa --te ALL  --chr Chr3 --number 200 --prefix MSU7.Chr3.ALL.rep3

echo "Simulate illumina reads"
python Simulation_illumina.py --input ./MSU7.Chr3.ALL --size 500

echo "FLY603_Chr2L_pogo"
echo "Simulate genome"
python Simulation_TE_Insertion.py --repeat ./Reference/flyte.fa --genome ./Reference/FLY603.Chr2L.fa --code FLY603 --te pogo  --chr Chr2L --number 200 --prefix FLY603.Chr2L.pogo.rep1
python Simulation_TE_Insertion.py --repeat ./Reference/flyte.fa --genome ./Reference/FLY603.Chr2L.fa --code FLY603 --te pogo  --chr Chr2L --number 200 --prefix FLY603.Chr2L.pogo.rep2
python Simulation_TE_Insertion.py --repeat ./Reference/flyte.fa --genome ./Reference/FLY603.Chr2L.fa --code FLY603 --te pogo  --chr Chr2L --number 200 --prefix FLY603.Chr2L.pogo.rep3

echo "Simulate illumina reads"
python Simulation_illumina.py --input FLY603.Chr2L.pogo --size 500

echo "TAIR10_Chr1_mPing"
echo "Simulate genome"
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/TAIR10.Chr1.fa --code TAIR10 --te mPing  --chr Chr1 --number 200 --prefix TAIR10.Chr1.mPing.rep1
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/TAIR10.Chr1.fa --code TAIR10 --te mPing  --chr Chr1 --number 200 --prefix TAIR10.Chr1.mPing.rep2
python Simulation_TE_Insertion.py --repeat ./Reference/mping.fa --genome ./Reference/TAIR10.Chr1.fa --code TAIR10 --te mPing  --chr Chr1 --number 200 --prefix TAIR10.Chr1.mPing.rep3

echo "Simulate illumina reads"
python Simulation_illumina.py --input TAIR10.Chr1.mPing --size 500


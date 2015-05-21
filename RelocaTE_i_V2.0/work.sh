#echo "MSU7.Chr4.mPing"
#python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa 
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/MSU7.Chr4.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
#python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing

#echo "MSU7.Chr4.ALL"
#python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/MSU7.Chr4.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
#python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL

#echo "MSU7.Chr3.mPing"
#python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/MSU7.Chr3.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
#python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.mPing

#echo "MSU7.Chr3.ALL"
python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/MSU7.Chr3.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL

#echo "FLY603.Chr2L.pogo"
#python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/FLY603.Chr2L.pogo --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/FLY603.Chr2L.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/pogo.fa
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/FLY603.Chr2L.pogo --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/FLY603.Chr2L.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/pogo.fa
#python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/FLY603.Chr2L.pogo

#echo "TAIR10.Chr1.mPing"
#python RunRelocaTEi.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/TAIR10.Chr1.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/TAIR10.Chr1.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
#python RunRelocaTEi_bam.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/ReadMapping/TAIR10.Chr1.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/TAIR10.Chr1.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa
#python SumReCall.py --call RelocaTEi --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/TAIR10.Chr1.mPing


#bash MSU7.Chr3.ALL.sh > MSU7.Chr3.ALL.sh.log 2>&1 &
#bash MSU7.Chr4.ALL.sh > MSU7.Chr4.ALL.sh.log 2>&1 &

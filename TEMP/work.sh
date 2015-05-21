echo "MSU7.Chr4.mPing"
python RunTEMP.py --input ../ReadMapping/MSU7.Chr4.mPing > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.mPing

echo "MSU7.Chr4.ALL"
python RunTEMP.py --input ../ReadMapping/MSU7.Chr4.ALL --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL

echo "MSU7.Chr4.ALL, call mPing only"
python RunTEMP.py --input ../ReadMapping/MSU7.Chr4.ALL

echo "MSU7.Chr3.mPing"
python RunTEMP.py --input ../ReadMapping/MSU7.Chr3.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping.fa > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.mPing

echo "MSU7.Chr3.ALL"
python RunTEMP.py --input ../ReadMapping/MSU7.Chr3.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL

echo "FLY603.Chr2L.pogo"
python RunTEMP.py --input ../ReadMapping/FLY603.Chr2L.pogo --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/FLY603.Chr2L.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping_pogo.fa > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/FLY603.Chr2L.pogo

echo "TAIR10.Chr1.mPing"
python RunTEMP.py --input ../ReadMapping/TAIR10.Chr1.mPing --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/TAIR10.Chr1.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/mping_pogo.fa > log 2> log2 &
python SumReCall.py --call TEMP --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/TAIR10.Chr1.mPing

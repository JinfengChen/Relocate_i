#echo "Mix.ALL.Repeat"
#python RunRelocaTE_multiTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Mix.ALL.Repeat --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/Mix.ALL.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/repeat.fa
#python SumReCall.py --call RelocaTE --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.mPing

#echo "MSU7.Chr4.ALL"
#python RunRelocaTE_multiTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
#python RunRelocaTE_multiTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr4.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.TSD.fa
#python SumReCall_multiTE.py --call RelocaTE --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr4.ALL

echo "MSU7.Chr3.ALL"
python RunRelocaTE_multiTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.fa
#python RunRelocaTE_multiTE.py --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa --repeat /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/RiceTE.TSD.fa
#python SumReCall_multiTE.py --call RelocaTE --input /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/MSU7.Chr3.ALL

################################################
# ADMIXTURE projection with harmonized dataset #
################################################

# Making a reference file based on P and K-DET bim
paste K-DET_1KGONL.bim K-DET_1KGONL.9.P | awk '{print $2" | "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15}'  | sort -k 1,1 > ReferenceFreqs.dat
awk '{print $2}' K-DET_1KGONL.bim | sort > ReferenceSNPsUsed.dat

# Subsetting and removing reference samples from earlier steps.
awk '{if ($6==30) print $1" "$2}' 3platform_1KGONL_PC_overlappingSNPs_RDY.fam > _3platform_NTR.dat

plink2 --bfile 3platform_1KGONL_PC_overlappingSNPs_RDY --extract ReferenceSNPsUsed.dat --keep _3platform_NTR.dat --make-bed --out 3platform_9POP

awk '{print $2" "$5" "$6" "NR}' 3platform_9POP.bim | sort -k 1,1 > _3platform_P1.dat
join -1 1 -2 1 _3platform_P1.dat ReferenceFreqs.dat | awk '{if ($2==$6 && $3==$7) print $4" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}'> _3platform_P2_A1A2.dat
join -1 1 -2 1 _3platform_P1.dat ReferenceFreqs.dat | awk '{if ($2==$7 && $3==$6) print $4" "(1-$8)" "(1-$9)" "(1-$10)" "(1-$11)" "(1-$12)" "(1-$13)" "(1-$14)" "(1-$15)" "(1-$16)}'> _3platform_P2_A2A1.dat
cat _3platform_P2_A1A2.dat _3platform_P2_A2A1.dat | sort -n -k 1,1 | awk '{print $2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10}' > 3platform_9POP.9.P.in

#CHECK
wc -l *_9POP.9.P.in
wc -l *_9POP.bim

rm $Admix/_*
rm $Admix/*.nosex

# Run ADMIXTURE projection
./admixture -j8 -P 3platform_9POP.bed 9
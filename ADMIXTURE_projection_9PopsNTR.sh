########################
# ADMIXTURE projection #
########################

# Making a reference file based on P and K-DET bim
paste K-DET_1KGONL.bim K-DET_1KGONL.9.P | awk '{print $2" | "$5" "$6" "$7" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15}'  | sort -k 1,1 > ReferenceFreqs.dat
awk '{print $2}' K-DET_1KGONL.bim | sort > ReferenceSNPsUsed.dat

# Subsetting and removing reference samples from earlier steps.
awk '{if ($6==30) print $1" "$2}' AFFY6_1KGONL_PC_SNPs_RDY.fam > _AFFY6_NTR.dat
awk '{if ($6==30) print $1" "$2}' AXIOM_1KGONL_PC_SNPs_RDY.fam > _AXIOM_NTR.dat
awk '{if ($6==30) print $1" "$2}' ILLGSA_1KGONL_PC_SNPs_RDY.fam > _ILLGSA_NTR.dat

plink2 --bfile AFFY6_1KGONL_PC_SNPs_RDY --extract ReferenceSNPsUsed.dat --keep _AFFY6_NTR.dat --make-bed --out AFFY6_9POP
plink2 --bfile AXIOM_1KGONL_PC_SNPs_RDY --extract ReferenceSNPsUsed.dat --keep _AXIOM_NTR.dat --make-bed --out AXIOM_9POP
plink2 --bfile ILLGSA_1KGONL_PC_SNPs_RDY --extract ReferenceSNPsUsed.dat --keep _ILLGSA_NTR.dat --make-bed --out ILLGSA_9POP

awk '{print $2" "$5" "$6" "NR}' AFFY6_9POP.bim | sort -k 1,1 > _AFFY6_P1.dat
join -1 1 -2 1 _AFFY6_P1.dat ReferenceFreqs.dat | awk '{if ($2==$6 && $3==$7) print $4" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}'> _AFFY6_P2_A1A2.dat
join -1 1 -2 1 _AFFY6_P1.dat ReferenceFreqs.dat | awk '{if ($2==$7 && $3==$6) print $4" "(1-$8)" "(1-$9)" "(1-$10)" "(1-$11)" "(1-$12)" "(1-$13)" "(1-$14)" "(1-$15)" "(1-$16)}'> _AFFY6_P2_A2A1.dat
cat _AFFY6_P2_A1A2.dat _AFFY6_P2_A2A1.dat | sort -n -k 1,1 | awk '{print $2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10}' > AFFY6_9POP.9.P.in

awk '{print $2" "$5" "$6" "NR}' AXIOM_9POP.bim | sort -k 1,1 > _AXIOM_P1.dat
join -1 1 -2 1 _AXIOM_P1.dat ReferenceFreqs.dat | awk '{if ($2==$6 && $3==$7) print $4" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}'> _AXIOM_P2_A1A2.dat
join -1 1 -2 1 _AXIOM_P1.dat ReferenceFreqs.dat | awk '{if ($2==$7 && $3==$6) print $4" "(1-$8)" "(1-$9)" "(1-$10)" "(1-$11)" "(1-$12)" "(1-$13)" "(1-$14)" "(1-$15)" "(1-$16)}'> _AXIOM_P2_A2A1.dat
cat _AXIOM_P2_A1A2.dat _AXIOM_P2_A2A1.dat | sort -n -k 1,1 | awk '{print $2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10}' > AXIOM_9POP.9.P.in

awk '{print $2" "$5" "$6" "NR}' ILLGSA_9POP.bim | sort -k 1,1 > _ILLGSA_P1.dat
join -1 1 -2 1 _ILLGSA_P1.dat ReferenceFreqs.dat | awk '{if ($2==$6 && $3==$7) print $4" "$8" "$9" "$10" "$11" "$12" "$13" "$14" "$15" "$16}'> _ILLGSA_P2_A1A2.dat
join -1 1 -2 1 _ILLGSA_P1.dat ReferenceFreqs.dat | awk '{if ($2==$7 && $3==$6) print $4" "(1-$8)" "(1-$9)" "(1-$10)" "(1-$11)" "(1-$12)" "(1-$13)" "(1-$14)" "(1-$15)" "(1-$16)}'> _ILLGSA_P2_A2A1.dat
cat _ILLGSA_P2_A1A2.dat _ILLGSA_P2_A2A1.dat | sort -n -k 1,1 | awk '{print $2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10}' > ILLGSA_9POP.9.P.in

#CHECK
wc -l *_9POP.9.P.in
wc -l *_9POP.bim

rm _*
rm *.nosex

# Run ADMIXTURE projection
admixture -j8 -P AXIOM_9POP.bed 9
admixture -j8 -P AFFY6_9POP.bed 9
admixture -j8 -P ILLGSA_9POP.bed 9


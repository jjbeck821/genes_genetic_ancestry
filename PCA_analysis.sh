################
# PCA analysis #
################

# INPUT DATA VARIABLES
AFFY6=AFFY6_CLN_SNP
AXIOM=AXIOM_CLN_SNP
ILLGSA=ILLGSA_CLN_SNP
KG1=1KG_P3V5A_B37
GONL=GONLV4

# Filtering data for PCAs.
king -b $AXIOM".bed" --unrelated --degree 2 --prefix AXIOM_
king -b $AFFY6".bed" --unrelated --degree 2 --prefix AFFY6_
king -b $ILLGSA".bed" --unrelated --degree 2 --prefix ILLGSA_

plink2 --bfile $AXIOM --keep AXIOM_unrelated.txt --make-founders --maf 0.01 --chr 1-22 --hwe 0.001 --geno 0.05 --mind 0.05 --make-bed --out AXIOM_UNR
plink2 --bfile $AFFY6 --keep AFFY6_unrelated.txt --make-founders --maf 0.01 --chr 1-22 --hwe 0.001 --geno 0.05 --mind 0.05 --make-bed --out AFFY6_UNR
plink2 --bfile $ILLGSA --keep ILLGSA_unrelated.txt --make-founders --maf 0.01 --chr 1-22 --hwe 0.001 --geno 0.05 --mind 0.05 --make-bed --out ILLGSA_UNR

plink2 --bfile AXIOM_UNR --indep-pairwise 250kb 1 0.50 --out AXIOM_PRU
plink2 --bfile AFFY6_UNR --indep-pairwise 250kb 1 0.50 --out AFFY6_PRU
plink2 --bfile ILLGSA_UNR --indep-pairwise 250kb 1 0.50 --out ILLGSA_PRU

plink2 --bfile AXIOM_UNR --extract AXIOM_PRU.prune.in --make-bed --out _TAX
plink2 --bfile AFFY6_UNR --extract AFFY6_PRU.prune.in --make-bed --out _TAF
plink2 --bfile ILLGSA_UNR --extract ILLGSA_PRU.prune.in --make-bed --out _TIL

# Notice plink1 here
plink1 --bfile _TAX --exclude long_range_LD_plink_format_build_37.txt --range --make-bed --out AXIOM_PC_SNPs --noweb
plink1 --bfile _TAF --exclude long_range_LD_plink_format_build_37.txt --range --make-bed --out AFFY6_PC_SNPs --noweb
plink1 --bfile _TIL --exclude long_range_LD_plink_format_build_37.txt --range --make-bed --out ILLGSA_PC_SNPs --noweb
rm _T*.bed
rm _T*.bim
rm _T*.fam


# Aligning this data with 1000G
awk '{print $2"_"$5"_"$6}' *_PC_SNPs.bim > _NA1A2.dat
awk '{print $2"_"$6"_"$5}' *_PC_SNPs.bim > _NA2A1.dat
cat _NA1A2.dat _NA2A1.dat | sort -u > _NPOSALLELES.dat
awk '{print $1"_"$4"_"$5"_"$6"\t"$2}' $KG1".bim" | sort -k1,1 > _N1KGALL.dat
join -1 1 -2 1 _NPOSALLELES.dat _N1KGALL.dat | sed 's/_/\t/g' | awk '{print $1"_"$2" "$5}' | sort -u > _INboth.dat
awk '{print $1}' _INboth.dat > IN1KG_SNPs_CHRPOS.dat
awk '{print $2}' _INboth.dat > IN1KG_SNPs_RS.dat
awk '{print $1" "$2}' _INboth.dat > IN1KG_SNPs_IDupdate.dat
rm _*.dat


# Aligning with GONL & 1KG and picking all unrelated in reference panels.
# Note since being unrelated is dependent on the input SNPs select the HM3 SNP subset.
awk '{print $5}' HM3_B37_map.dat > _HM3_RSs.dat
awk '{print $1"_"$2"_"$3"_"$4}' HM3_B37_map.dat > _HM3_4POS.dat
awk '{print $1"_"$2"_"$4"_"$3}' HM3_B37_map.dat >> _HM3_4POS.dat
plink2 --bfile $GONL --extract _HM3_4POS.dat --make-bed --out _GONLHM3
plink2 --bfile $KG1 --extract _HM3_RSs.dat --make-bed --out _1KGHM3
king -b _GONLHM3.bed --unrelated --degree 2 --prefix GONL_
king -b _1KGHM3.bed --unrelated --degree 2 --prefix 1KG_
plink2 --bfile $KG1 --extract IN1KG_SNPs_RS.dat --make-bed --out 1KG_PC_SNPs_RS
awk '{print $1"_"}' IN1KG_SNPs_CHRPOS.dat > _Find.dat
awk '{print $2}' $GONL".bim" | grep -f _Find.dat > INGONL_SNPs_start.dat
rm _Find.dat
plink2 --bfile $GONL --extract INGONL_SNPs_start.dat --make-bed --out GONL_PC_SNPs_CHRBPA1A2
plink2 --bfile GONL_PC_SNPs_CHRBPA1A2 --keep GONL_unrelated.txt --make-founders --make-bed --out GONL_PC_SNPs_UNR_CHRBPA1A2
mv GONL_PC_SNPs_UNR_CHRBPA1A2.bim GONL_PC_SNPs_UNR_CHRBPA1A2.bim_org
awk '{print $1" "$1"_"$4" "$3" "$4" "$5" "$6}' GONL_PC_SNPs_UNR_CHRBPA1A2.bim_org > GONL_PC_SNPs_UNR_CHRBPA1A2.bim
plink2 --bfile GONL_PC_SNPs_UNR_CHRBPA1A2 --update-name IN1KG_SNPs_IDupdate.dat --make-bed --out GONL_PC_SNPs_UNR_RS
awk '{print $2}' GONL_PC_SNPs_UNR_RS.bim | grep -v "rs" > GONL_remove_these_snps_still.dat
plink2 --bfile GONL_PC_SNPs_UNR_RS --exclude GONL_remove_these_snps_still.dat --make-bed --out GONL_PC_SNPs_UNR_RS_1
cat GONL_PC_SNPs_UNR_RS_1.bim 1KG_PC_SNPs_RS.bim | awk '{print $2}' | sort | uniq -c | awk '{if ($1==2) print $2}' > INBothRefsRS.dat
plink2 --bfile 1KG_PC_SNPs_RS --keep 1KG_unrelated.txt --make-bed --out 1KG_PC_SNPs_UNR_RS
plink2 --bfile 1KG_PC_SNPs_UNR_RS --bmerge GONL_PC_SNPs_UNR_RS_1 --extract INBothRefsRS.dat --geno 0.02 --make-bed --out ./Saved/1KGnGONL_PC_SNPs_REF

# Generating PCA datasets
awk '{print $2}' AXIOM_PC_SNPs.bim > AXIOM_SNP_List.dat
awk '{print $2}' AFFY6_PC_SNPs.bim > AFFY6_SNP_List.dat
awk '{print $2}' ILLGSA_PC_SNPs.bim > ILLGSA_SNP_List.dat

plink2 --bfile $AXIOM --extract AXIOM_SNP_List.dat --make-bed --out AXIOM_FullSample_PC_SNPs_1 
plink2 --bfile $AFFY6 --extract AFFY6_SNP_List.dat --make-bed --out AFFY6_FullSample_PC_SNPs_1
plink2 --bfile $ILLGSA --extract ILLGSA_SNP_List.dat --make-bed --out ILLGSA_FullSample_PC_SNPs_1

plink2 --bfile AXIOM_FullSample_PC_SNPs_1 --update-name IN1KG_SNPs_IDupdate.dat --extract INBothRefsRS.dat --make-bed --out AXIOM_FullSample_PC_SNPs_2
plink2 --bfile AFFY6_FullSample_PC_SNPs_1 --update-name IN1KG_SNPs_IDupdate.dat --extract INBothRefsRS.dat --make-bed --out AFFY6_FullSample_PC_SNPs_2
plink2 --bfile ILLGSA_FullSample_PC_SNPs_1 --update-name IN1KG_SNPs_IDupdate.dat --extract INBothRefsRS.dat --make-bed --out ILLGSA_FullSample_PC_SNPs_2

plink2 --bfile AXIOM_FullSample_PC_SNPs_2 --bmerge ./Saved/1KGnGONL_PC_SNPs_REF --geno 0.02 --make-bed --out AXIOM_1KGONL_PC_SNPs_1
plink2 --bfile AFFY6_FullSample_PC_SNPs_2 --bmerge ./Saved/1KGnGONL_PC_SNPs_REF --geno 0.02 --make-bed --out AFFY6_1KGONL_PC_SNPs_1
plink2 --bfile ILLGSA_FullSample_PC_SNPs_2 --bmerge ./Saved/1KGnGONL_PC_SNPs_REF --geno 0.02 --make-bed --out ILLGSA_1KGONL_PC_SNPs_1

# Cleaning up and avoiding GONL overlap adding expected ancestry
king -b AXIOM_1KGONL_PC_SNPs_1.bed --duplicate --prefix AX_Dupes
king -b AFFY6_1KGONL_PC_SNPs_1.bed --duplicate --prefix AF_Dupes
king -b ILLGSA_1KGONL_PC_SNPs_1.bed --duplicate --prefix IG_Dupes

grep "A" *.con | sed 's/:/\t/' | awk '{print $2}' | sort -u > FAMinGONL.dat
grep -w -f FAMinGONL.dat *_PC_SNPs_1.fam | sed 's/:/ /' | awk '{print $2" "$3}' > GONL_out.dat

cat GONL_out.dat Others_no_consent_out.dat | sort -u > All_out.dat

cat 1KG_ances_pheno.dat DATA_ances_pheno.dat GONL_ances_pheno.dat > ALL_ances_pheno.dat

plink2 --bfile AXIOM_1KGONL_PC_SNPs_1 --remove All_out.dat --pheno ALL_ances_pheno.dat --make-bed --out AXIOM_1KGONL_PC_SNPs_RDY --allow-no-sex
plink2 --bfile AFFY6_1KGONL_PC_SNPs_1 --remove All_out.dat --pheno ALL_ances_pheno.dat --make-bed --out AFFY6_1KGONL_PC_SNPs_RDY --allow-no-sex
plink2 --bfile ILLGSA_1KGONL_PC_SNPs_1 --remove All_out.dat --pheno ALL_ances_pheno.dat --make-bed --out ILLGSA_1KGONL_PC_SNPs_RDY --allow-no-sex

plink2 --bfile AXIOM_1KGONL_PC_SNPs_RDY --recode --out AXIOM_1KGONL_PC_SNPs_RDY --allow-no-sex
plink2 --bfile AFFY6_1KGONL_PC_SNPs_RDY --recode --out AFFY6_1KGONL_PC_SNPs_RDY --allow-no-sex
plink2 --bfile ILLGSA_1KGONL_PC_SNPs_RDY --recode --out ILLGSA_1KGONL_PC_SNPs_RDY --allow-no-sex

# Run PCA
smartpca -p AXIOMvs1KGONL.par
smartpca -p AFFY6vs1KGONL.par
smartpca -p ILLGSAvs1KGONL.par

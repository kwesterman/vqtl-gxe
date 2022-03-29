ANNOVARDIR=../opt/annovar
OUTPUTDIR=../data/processed/vqtl_ss/annovar
TARGETDIR=../data/processed/vqtl_ss/annovar

perl $ANNOVARDIR/annotate_variation.pl \
        -out $OUTPUTDIR/vqtl_variants \
        -build hg19 \
        $TARGETDIR/vqtl_variants.avinput \
        $ANNOVARDIR/humandb/

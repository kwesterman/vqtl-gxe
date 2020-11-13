ANNOVARDIR=../opt/annovar
OUTPUTDIR=../data/processed/vqtl_ss/annovar
TARGETDIR=../data/processed/vqtl_ss/annovar

perl $ANNOVARDIR/annotate_variation.pl \
        -out $OUTPUTDIR/vqtl_MA_hits \
        -build hg19 \
        $TARGETDIR/vqtl_MA_hits.avinput \
        $ANNOVARDIR/humandb/

rm tmp.avinput

#!/bin/sh

# association study 001
bin/plink_1.9_mac --file data/assoc_01 --allow-no-sex --assoc
bin/plink_1.9_mac --file data/assoc_01 --allow-no-sex --fisher
bin/plink_1.9_mac --file data/assoc_01 --allow-no-sex --model

cat plink.assoc
cat plink.assoc.fisher
cat plink.model

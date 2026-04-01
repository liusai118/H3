The Dataset folder is available for download at https://drive.google.com/file/d/1B_CFm9nzzaNJ7wFXYJAxSBadj3dUFNnh/view?usp=sharing
| 文稿中的句子/数字 | 脚本 | 输出文件 | 对应数值/说明 |
|------------------|------|----------|----------------|
| 11,048 positive sequences and 11,047 valid negative sequences | ./proprecess/gkmsvm/preprocess.R | ./proprecess/gkmsvm/all.csv | 22095 rows |
| model-specific thresholds | ./run/rep*.ipynb | ./out/threshold_manifest.json | rep1~rep10 |
| 3,000,000 candidate sequences | ./xai/random_02.R | ./xai/pos.csv | 3000000 rows |
| 100,000 mutated sequences | ./xai/mut_04.R | ./xai/mutant.csv | 100000 rows |
| 500000 sequences based on low-frequency 5-mers | ./xai/random_02.R | ./xai/neg.csv | 600000 rows |
| 10,000 optimized decoy conformations | ./rosetta/data_analysis_04.ipynb | ./Dataset/rosetta.csv | 10000 rows |
| 10 sequence | ./rosetta/rosetta_05.R | ./rosetta/result.csv | 10 rows |

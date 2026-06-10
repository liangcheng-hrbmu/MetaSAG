# MetaSAG Usage 
## Step 11. MetaK-Lytic

This is an algorithm for predicting the lysis ability of bacteriophages based on 31-mer sequence data

## 1. art-illumina 

Simulate the sequencing results of the second-generation sequencing in this way, just run

```bash
art_illumina -ss HS20 -i your_sequence.fasta -l 100 -f 20 -o your_output_file
```

## 2. k-mer sequence frequency feature extraction, just run

```bash
k-mer.py
```

Obtain the sequence features in the.npy format

## 3. Integrate all the feature frequency files, just run

```bash
k_mer_union.py
```

## 4. Make predictions using MetaK_Lytic

```bash
MetaK_Lytic.py --mode train --temperate-file temperate_phgaes_features.npy --lytic-file lytic_phages_features.npy

MetaK_Lytic.py --mode test --temperate-file ./temperate_features.npy --lytic-file ./lytic_features.npy
```

## 5. pre-trained model

If you don't want to use your own data for training, here we provide you with a pre-trained model named "final_model.h5".
Note that at this point, <mark>the first three steps still need to be followed to obtain the.npy format sequence feature file</mark>.


Download it, and then run `phage_predictor.py`


## Pre-trained Model
- [`final_model.h5`](../Example_data/final_model.h5)

## Example Files
1. [`example_Lytic.npy`](../Example_data/example_Lytic.npy)  
2. [`example_temperate_feature.npy`](../Example_data/example_temperate_feature.npy)
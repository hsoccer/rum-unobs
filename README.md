<!---
contributors:
  Haruki Kono
  Alec Sandroni
--->

# Replication Package for "Random Utility with Unobservable Alternatives"

This folder contains all code and data required to reproduce the results in Section V of “Random Utility with Unobservable Alternatives” by Haruki Kono, Kota Saito, and Alec Sandroni. In particular, the code computes the lower and upper bounds on unobservable choice probabilities implied by the random utility model using the dataset provided by McCausland et al. (2020). We also include code to generate the aggregated dataset from the McCausland et al. (2020) data.

## 1. Data Availability Statement

The analysis executed in this folder is based on the dataset obtained in an experimental study by McCausland et al. (2020).
This is available at ./Data_Tables/RCM_multi_data.tex in their supplementary material (https://doi.org/10.1093/ej/uez039).
This dataset is provided in the form of .tex. 


## 2. Computational Requirements

- macOS Monterey 12.3.1
- Python 3.9.7
    - numpy 1.26.4
    - pandas 1.5.3
    - matplotlib 3.7.1
    - scipy 1.7.1
    - the file ./requirements.txt lists these dependencies, please run `pip install -r requirements.txt` as the first step.
- Expected running time is less than 10 seconds.


## 3. Overview of Directory Structure

### . (Top directory)

- README: this file
- collective_choice_data.csv: a csv containing the choice frequencies calculated from the raw McCausland et al. (2020) dataset.
- make_lottery_figure.py: a code that makes the figure named ./id_comparison.png that appears in Section V of the paper.
- requirements.txt

### ./packages

- bm.py: a set of functions related to the computation of Block-Marschak polynomials
- utils.py: a set of other useful functions


## 4. Steps to Reproduce Results in Paper

## 4.1 Steps to reproduce the figure

To create the figure from Section V of the paper, execute ./make_lottery_figure.py, which can be done in Unix by issuing the command

```python
python make_lottery_figure.py
```
This script (i) constructs the calibrated dataset from the choice frequencies in collective_choice_data.csv, (ii) computes the implied lower and upper bounds on unobservable choice probabilities under the random utility model, and (iii) generates the figure reported in Section V of the paper.

Section 4.2 describes how to generate the aggregated dataset collective_choice_data.csv. Reproducing collective_choice_data.csv is not required to reproduce the figure, because the aggregated dataset is included in this replication package, courtesy of the authors of McCausland et al. (2020).

## 4.2 Steps to reproduce aggregated dataset (collective_choice_data.csv)

To reproduce the aggregated dataset, first get RCM_multi_data.tex. This is available at ./Data_Tables/RCM_multi_data.tex in the supplementary material of McCausland et al. (2020) (https://doi.org/10.1093/ej/uez039). Make sure that RCM_multi_data.tex is in the folder contaning the code.

Then execute ./from_tex_to_csv.py, which can be done in Unix by issuing the command

```python
python from_tex_to_csv.py
```

This writes the processed individual data to the Raw Data folder. Next, execute ./process_data.py by issuing the command

```python
python process_data.py
```

This creates the collective_choice_data.csv file, which contains the aggregated choice frequencies for each choice sets.

## 5. Results

After executing ./make_lottery_figure.py, ./outputs/id_comparison.png will be saved.
The image should look like:

![id_comparison](./outputs/id_comparison.png)

## 6. References

- McCausland, W. J., C. Davis-Stober, A. A. Marley, S. Park, and N. Brown (2020): "Testing the random utility hypothesis directly," The Economic Journal, 130(625), 183–207.

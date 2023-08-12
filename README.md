---
contributors:
  - Haruki Kono
  - Kota Saito
  - Alec Sandroni
---

# REPLICATION PACKAGE FOR "Axiomatization of Random Utility Model with Unobservable Alternatives"

This folder provides all the codes and data to reproduce the results Section 4 of the paper "Axiomatization of Random Utility Model with Unobservable Alternatives" by Haruki Kono, Kota Saito, Alec Sandroni (https://arxiv.org/abs/2302.03913).
In particular, the codes compute the lower and upper bounds of unobservable choice probabilities implied by the random utility model using the dataset provided by McCausland et al. (2020).

## 1. Data Availability Statement

The analysis executed in this folder is based on the dataset obtained in an experimental study by McCausland et al. (2020).
This is available at ./Data_Tables/data_tables.pdf in their supplementary material (https://doi.org/10.1093/ej/uez039).
This dataset is provided in the form of .pdf (and .tex), but we translated it into .csv files in ./Raw Data directory.


## 2. Computational Requirements

- macOS Monterey 12.3.1
- Python 3.9.7
    - numpy 1.24.3
    - pandas 1.3.4
    - matplotlib 3.7.1
    - scipy 1.7.1
    - the file ./requirements.txt lists these dependencies, please run `pip install -r requirements.txt` as the first step.
- Expected running time is less than 10 seconds.


## 3. Overview of Directory Structure

### . (Top directory)

- README: this file
- process_data.py: a code that aggregates all files in Raw Data directory, create a single stochastic choice dataset, and save it as ./collective_choice_data.csv
- make_lottery_figure.py: a code that makes the figure named ./id_comparison.png that appears in Section 4 of the paper


### ./Raw Data

This directory contains the results of the experimental study by McCausland et al. (2020).
In the experiment, the authors prepared five different lotteries and asked 141 subjects to choose a lottery from each subset of the five.
They repeated this six times for each individual.

This directory contains 141 files and each of them has each participant's results.
The numbers in the column (0 - 4) represent the five lotteries. 
The numbers in the row (0 - 25) represent all possible choice sets in the following correspondance:

| Index |     Choice set     |
|:-----:|:------------------:|
|    0  |        (0, 1)      |
|    1  |        (0, 2)      |
|    2  |        (1, 2)      |
|    3  |      (0, 1, 2)     |
|    4  |        (0, 3)      |
|    5  |        (1, 3)      |
|    6  |      (0, 1, 3)     |
|    7  |        (2, 3)      |
|    8  |      (0, 2, 3)     |
|    9  |      (1, 2, 3)     |
|   10  |     (0, 1, 2, 3)   |
|   11  |        (0, 4)      |
|   12  |        (1, 4)      |
|   13  |      (0, 1, 4)     |
|   14  |        (2, 4)      |
|   15  |      (0, 2, 4)     |
|   16  |      (1, 2, 4)     |
|   17  |     (0, 1, 2, 4)   |
|   18  |        (3, 4)      |
|   19  |      (0, 3, 4)     |
|   20  |      (1, 3, 4)     |
|   21  |     (0, 1, 3, 4)   |
|   22  |      (2, 3, 4)     |
|   23  |     (0, 2, 3, 4)   |
|   24  |     (1, 2, 3, 4)   |
|   25  |   (0, 1, 2, 3, 4)  |

The cells in the table represent the frequencies of choice.
If a lottery is unavailable at a choice set, then the corresponding cell is left blank.


### ./packages

- bm.py: a set of functions related to the computation of Block-Marschak polynomials
- utils.py: a set of other useful functions


## 4. Steps to Reproduce Results in Paper

It is assumed that the working directory is ./.

(a) Process the raw dara.

First, aggregate 141 files in ./Raw Data to create ./collective_choice_data.csv.
To do so, execute ./process_data.py, which can be done in Unix by issuing the command

> python process_data.py


(b) Make the figure.

Next, compute the lower and upper bounds of unobservable choice probabilities implied by the random utility model and draw the figure in Section 4 of the paper.
To do so, execute ./make_lottery_figure.py, which can be done in Unix by issuing the command

> python make_lottery_figure


## 5. Results

After executing ./make_lottery_figure.py, ./id_comparison.png will be saved.
The image should look like:

![id_comparison](./id_comparison.png)

## 6. References

- McCausland, W. J., C. Davis-Stober, A. A. Marley, S. Park, and N. Brown (2020): "Testing the random utility hypothesis directly," The Economic Journal, 130(625), 183–207.
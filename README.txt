The repository (https://github.com/hnasu-code/multi_wiener) contains the following files:

1. "run_Table5.R": The main R file to reproduce Table 5 in the paper.

2. "func_start.R": The R function that finds the helpful starting point for the EM algorithm.

3. "func_M0.R": The R function that performs the EM algorithm for model M0.

4. "func_M2.R": The R function that performs the EM algorithm for model M2.

5. "func_M3.R": The R function that performs the EM algorithm for model M3.

6. "func_M4.R": The R function that performs the EM algorithm for model M4.

7. "IRLED.csv": The CSV file storing the raw IRLED degradation Data.

Instructions: To execute the code, ensure that all the files mentioned above are placed in a common folder. Then, open "run_Table5.R" in either RStudio or R terminal and install all the relevant packages as claimed at the beginning. Run the code line by line under each code section in terms of a corresponding candidate model to obtain the results, which include the log-likelihood and AIC values.

# Sample Code for M-LCBBO

## Paper Source: "[Constrained Bayesian Optimization with Lower Confidence Bound](https://doi.org/10.1080/00401706.2024.2336535)"

### Instructions to run the code together.

1. <b>Requirement:</b>
    - R-version:>=4.1.3
        Dependency packages: tgp, laGP
    - Python-version:>=3.9
        Dependency package: scikit-learn:1.4.1.post1, numpy, scipy
 
    <t>Note: Before proceeding further, ensure you have installed all necessary dependencies in a single environment (R and Python).

2. <b>Running the R Script in the terminal:</b>

    - Navigate the folder containing the `figure3.R` and `figure3.py` files.
    
    - Execute the R file using `Rscript figure3.R`.<br>
            <t>Note: This will run the `figure3.py` file as well (internally).

    - Provide the values for number of initial samples, number of additional samples, and number of Monte Carlo iterations as arguments:<br>
    
    ...Default values: (a) Number of initial samples: 20<br>
                            (b) Number of additional samples: 100<br>
                            (c) Number of Monte Carlo runs: 30<br>
            Note: a. To run with default values simply type "0" for all the arguments, or else Provide integer values as per your requirement.
                  b. (i). Number of initial samples should be greater than or equal to (>=) 10.
                     (ii). Number of additional samples must be greater than number of initial samples.
                     (iii). Number of Monte Carlo runs is independent of other two parameters.

    IV. Wait to finish all computations (This may take some time according to your arguments) to see the Plot.jpeg.


3. <b>If the above step 2 doesn't work then try the following steps:<b>

    I. Run the python file first using the command `python figure3.py` of "python3 figure3.py" with arguments.
            ## Example:  "python figure3.py 10 30 3" in the terminal. (Arguments are separated by space.)

    II. Open "figure3.R" in an R IDE(like RStudio), and comment out lines 345 to 358, where Python script is called using system command.
    
    III. Now Execute the R file using the command "Rscript figure3.R" and then provide same arguments (given in Python).
            ## Example:  `Rscript figure3.R` in the terminal. 
                        -- This will provoked step 2(III) in command line.
    
    IV. Wait for some time until you see the Plot.jpeg.

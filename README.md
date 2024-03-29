# Sample Code for M-LCBBO

### This is a sample code for the article: "[Constrained Bayesian Optimization with Lower Confidence Bound](https://doi.org/10.1080/00401706.2024.2336535)"

### Instructions to run the code together.

1. <b>Requirement:</b>
    1. R-version: $\geq$ 4.1.3
        Dependency packages: tgp, laGP
    2. Python-version: $\geq$ 3.9
        Dependency package: scikit-learn:1.4.1.post1, numpy, scipy
 
    <t>Note: Before proceeding further, ensure you have installed all necessary dependencies in a single environment (R and Python).

2. <b>Running the R Script in the terminal:</b>

   1. Navigate the folder containing the `figure3.R` and `figure3.py` files.

   2. Execute the R file using `Rscript figure3.R`.<br>

      <t>Note: This will run the `figure3.py` file as well (internally).
   
   3. <a id="third-heading"> </a> Provide the values for number of initial samples, number of additional samples, and number of Monte Carlo iterations as arguments:<br>
        - Default values:<br>
            - Number of initial samples: 20
            - Number of additional samples: 100
            - Number of Monte Carlo runs: 30
        - <t>Notes: <br>
            - To run with default values simply type "0" for all the arguments, or else Provide integer values as per your requirement.<br>
            - Number of initial samples should be greater than or equal to ($\geq$) 10.<br>
            - Number of additional samples must be greater than number of initial samples.<br>
            - Number of Monte Carlo runs is independent of other two parameters.<br>
    
    4. Wait to finish all computations (This may take some time according to your arguments) to see the `Plot.jpeg`.


4. <b>If the above step 2 doesn't work then try the following steps:</b>

    1. Run the Python file first using the command `python figure3.py` or `python3 figure3.py` with arguments.<br>
        - Example:  `python figure3.py 10 30 3` in the terminal. (Arguments are separated by space.)

    2. Open `figure3.R` in an R IDE (like RStudio), and comment out lines 345 to 358, where Python script is called using system command.
    
    3. Now Execute the R file using the command `Rscript figure3.R` and then provide same arguments (given in Python).
        - Example:  `Rscript figure3.R` in the terminal.
            - This will provoke [step 2(iii)](#third-heading) in command line.
    
    4. Wait for some time until you see the Plot.jpeg.

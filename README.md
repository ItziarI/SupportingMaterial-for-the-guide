# Association Analysis Between Polygenic Risk Scores and Traits: Practical Guidelines and Tutorial with an Illustrative Data Set of Schizophrenia

Supported Information of the tutorial paper *Association Analysis Between Polygenic Risk Scores and Traits: Practical Guidelines and Tutorial with an Illustrative Data Set of Schizophrenia* can be reproduced with the provided material in the repository. 

Examples of different complexity are also included to better understand the underlying steps and procedure. In all cases, the raw data and the R code are available so that the reader can reproduce its own output. Nevertheless, the generated *pdf* files have also been included to ease the reading, if wanted.


 
**Included examples with given data and R codes**

   - Sinthetics: 
   
       - *WExample1.csv* and *WorkingExample1_code.Rmd*: 
          Example with continuous trait, and model fulfilling linearity, normality and constant variance assumptions.
       
       - *WExample2.csv* and *WorkingExample2_code.Rmd*: Example with continuous trait and the steps taken to address issues derived from the assumption of normality of errors.
       
       - *WExample3.csv* and *WorkingExample3_code.Rmd*: Example with continuous trait and the steps taken to address issues in the initial fitted model with non-constant variance.
       
       - *WExample4.csv* and *WorkingExample4_code.Rmd*: Example with binary trait. Therein, a factor, with different association strengths for its levels is handled. 
   

   - Real data set: 
       + *Real_data_Negative.csv* with codes in *Real_dataCAPE_Negative_code.Rmd* 
       + *Real_data_Positive.csv* with codes in  *Real_dataCAPE_Positive_code.Rmd*
      
           Data set contains PRSs for psychotic-like experiences (PLE) measured on 227 healthy individuals. PLE are similar psychotic experiences to those experienced by patients with schizophrenia but found in an attenuated form in healthy subjects. PLE are considered to be normally distributed in the general population, with just a few individuals presenting high levels of PLE and thus being the ones at risk of developing psychosis. 
   

   

 
 **To run the examples**
 
   Once all files are download, each .Rmd file can be executed/knitted in R and a *pdf* file with the corresponding results is generated.
   
   Nevertheless, for the reader that do not want to execute the code and to ease the reading of the tutorial, the final pdf files are also included.


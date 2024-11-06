The repository contains programs to obtain the expression of Lyapunov density as per Theorem 1 in the paper "Construction of Lyapunov Certificates for Systems on Torus Using Trigonometric Polynomials" (https://www.preprints.org/manuscript/202409.1361/v1).

Step 1- Use the Mathematica file "Step1-Compute_Fourier_Coefficients_of_Oscillator_System.nb". The changes that might be reflected with different examples are:
                  i... the number and expressions for component functions
                  ii.. the number of harmonics
                  iii. the number of phase (or phase-difference) variables in the component functions
        This step outputs a matrix where each row consists of an index followed by the coefficient of the component function corresponding to the index. 
        
Step 2- Input the matrix obtained in the previous case to the MATLAB file "F.m". For computational efficiency, you may skip. 
                  i.. the rows where the second half of the entries are zero. For instance, you may skip the row (1 2 3 0 0 0) while inputting in "F.m".
                  ii. the rows whose first half is negative of the first half of a previous row. For instance, suppose the matrix has the following rows 
                                                                      ( *  *  * * * *) : * meaning any arbitrary element
                                                                      (-1 -2 -3 a b c)
                                                                      ( 1  2  3 x y z)
                                                                      ( *  *  * * * *)
                      Then, you may choose to skip typing the third row since the first half of the third row is -1 times the first half of the second row. 
         If the number of harmonics increases or decreases, corresponding changes should be reflected in MATLAB files (refer to individual files for information on achieving that). Run "driver.m" after making all            the changes.

Step 3- When the MATLAB code returns "solved" as an output, save the Gram Matrix as "MatrixGV.mat" and import it to the Mathematica file "Step3-MATLAB_Output_to_Lyapunov_Function.nb," which returns the trigonometric polynomial V such that 1/V is a Lyapunov density. Don't forget to make changes in the Mathematica file on the change of dimensions (refer to the file for details).
                                                                                                                                                                

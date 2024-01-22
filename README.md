This Jupyter Notebook and Python files will create an analytic sounding that can be used for submission to run CM1. 
An Anaconda environment will need to be set up with the Metpy and Sharppy modules installed. 

There are quite a few comments throughout the code to help users understand what each portion is doing. 

This current edition has the thermodynamic profile generator thanks to Robert Warren, who provided his code for creating analytic soundings in his paper "Impact of Variations in Upper-Level Shear on 
Simulated Supercells" Warren et al. 2017. This analytic profile code is also largely based off of the McCaul and Weisman 2001. All code provided by Robert was in IDL - this is just converted to Python. 
Therefore, I take zero credit for myself :) 

This current edition for the kinematic profile generator is fairly crude and developed by myself. I know there are a lot of analytic hodograph generating codes out there, therefore I would only recommend using
this if you need a very simple hodograph, as it works just sufficiently for a straight hodograph or a quarter-circle hodograph. 

This code also has the option to save the figure of the profile along with the text file that is saved for CM1 submission.

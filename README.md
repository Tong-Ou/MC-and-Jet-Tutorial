# MC-and-Jet-Tutorial
This repository is committed to share the codes written by Tong (in Python) for the practice in Monte Carlo and Jet Tutorial.
All of the codes (except the plot_and_fit.py) are independent, so you can run any one of them immediately.
About the display of results: Except PartonShowers_2partons.py and JetClustering_incl/excl.py, the results of which are stored in txt files for later analysis, results of other programs are shown by histograms (some with fitting).
There are three PartonShowers_xx.py in the repository. They are designed for different tasks: (a) The PartonShowerModel.py is the simplest one and is designed to generate a parton shower in 2D space. (b) The PartonShowerModel_re.py is also designed to generate parton shower in 2D space but with a different PDF(\theta). (c) The PartonShowers_2partons.py can generate parton showers for two initial partons in 3D space.
The plot_and_fit.py is designed to analyze results of JetClustering_incl/excl.py. Of course you can run it for another set of data with some necessary modifications of the codes.
Last but not least, all of the codes are completed by Tong independently, misunderstanding of physics process or improper design of codes may occur in the programs. Please let me know whenever you find some bugs or have some suggestions for improvement. Thanks for your patience to a novice :)

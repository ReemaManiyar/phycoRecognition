# phycoRecognition

Dependencies:
1. Python 3.6.0
2. Anaconda 4.3.1 64-bit
3. BioPython 1.67

To use phycoRecognition tool:
* Extract features:
    * Use src/feature.py file to extract features
    * This script will internally call various scripts from 'src' directory
* Load trained model
    * Use model/phycoTool.pkl file using joblib.load. 
    * This will return an object to sklearn.ensemble.RandomForestClassifier class
    * This object can be used to various class methods.
    * For reference: http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
* Run the test samples
    * Use 'predict' function to get the class of test contig
    * Alternatively use 'score' function to get the prediction probability
 

## Synopsis
[DNA shuffling](https://en.wikipedia.org/wiki/DNA_shuffling) is a molecular biology technique used in directed evolution studies.  The technique relies on the similarities between two or more DNA sequences.  This program can be used to increase the similarities between two DNA sequences.  After harmonization the two sequence can be shuffled together and screening for improved properties.    

This code was inspired by the paper [Shuffle Optimizer: A Program to Optimize DNA Shuffling for Protein Engineering](https://link.springer.com/protocol/10.1007%2F978-1-4939-6343-0_3) written by John N. Milligan and Daniel J. Garry.  However, the code below uses Python 3, and makes a few different codon choices.  On the whole though the outputs should be largely the same.


## Code Example
Run code with: 

```
 
# print sequences into the command prompt  
python shuffle.py example_sequences.fna 

# or wirte the outputs to a file
python shuffle.py example_sequences.fna > output.txt
```

## Installation

This program uses the Biopython library. Go [here](http://biopython.org/wiki/Download) for installation instructions before using the script.

## License

MIT License

Copyright (c) \[2017\] \[Jeremy LaBarge\]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
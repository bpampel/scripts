This contains some python scripts to calculate the wavelet coefficients and functions and do som initial tests.

The main work is the *"wavelet_generator.py"*, that can be used to generate coefficients or the form of the wavelets.
While this is fine to generate Db type wavelets and the resulting coefficients are in the Plumed code, this does not result in the canonical choice for the Sym type wavelets.

Therefore the *"sym_coefficients_permutations.py"* script was written to check if the literature values can actually be produced with my code: Yes they can. But the choice of the "right" coefficients for "most linear phase" appears to be not trivial, so I rather took the literature values for Plumed.

The "determine_wavelet_tails.py"* script calculates how much of the wavelets can be cut off to keep the range with all values that have at least 1% of the maximum value.
This is just for historic reasons and was implemented differently in Plumed.

# hydrogenic_matrix_elements
This is a python script that computes matrix elements of the form:

<img src="https://latex.codecogs.com/svg.latex?\int&space;d&space;\Vec{r}&space;d&space;\Vec{r'}&space;\&space;\frac{&space;\&space;\psi^*_{n_a,l_a,m_a}&space;(\Vec{r})&space;\&space;\psi^*_{n_c,l_c,m_c}&space;(\Vec{r'})&space;\&space;\psi_{n_b,l_b,m_b}&space;(\Vec{r})&space;\&space;\psi_{n_d,l_d,m_d}&space;(\Vec{r'})&space;}{|\Vec{r}&space;-&space;\Vec{r'}|}" title="\int d \Vec{r} d \Vec{r'} \ \frac{ \ \psi^*_{n_a,l_a,m_a} (\Vec{r}) \ \psi^*_{n_c,l_c,m_c} (\Vec{r'}) \ \psi_{n_b,l_b,m_b} (\Vec{r}) \ \psi_{n_d,l_d,m_d} (\Vec{r'}) }{|\Vec{r} - \Vec{r'}|}" />

Where <img src="https://latex.codecogs.com/svg.latex?\psi_{n,l,m}(\Vec{r})" title="\psi_{n,l,m}(\Vec{r})" /> is a hydrogen-like atomic orbital with charge Z. These matrix elements occur, for example, in Hartree-Fock theory for a multi electron atom when a hydrogen-orbital basis is used.

The program is designed to be as simple as possible to use. For example, to generate all possible matrix elements up to principal quantum number n=2 and write them to a file called 'outfile.csv', type:

python3 hydrogenic_matrix_elements.py -n_cutoff 2 -write_to_file outfile

Examples can be found in the 'examples.txt' file.

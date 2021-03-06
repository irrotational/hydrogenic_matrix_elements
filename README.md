# hydrogenic_matrix_elements
This is a python script that computes matrix elements of the form:

<img src="https://latex.codecogs.com/svg.latex?\bigg<&space;\&space;a(\Vec{r})&space;c(\Vec{r'})&space;\&space;\bigg|&space;\&space;\frac{1}{\Vec{r}-\Vec{r'}}&space;\&space;\bigg|&space;\&space;b(\Vec{r})&space;d(\Vec{r'})&space;\&space;\bigg>&space;\equiv&space;\int&space;d&space;\Vec{r}&space;d&space;\Vec{r'}&space;\&space;\frac{&space;\&space;\psi^*_{n_a,l_a,m_a}&space;(\Vec{r})&space;\&space;\psi^*_{n_c,l_c,m_c}&space;(\Vec{r'})&space;\&space;\psi_{n_b,l_b,m_b}&space;(\Vec{r})&space;\&space;\psi_{n_d,l_d,m_d}&space;(\Vec{r'})&space;}{|\Vec{r}&space;-&space;\Vec{r'}|}" title="\bigg< \ a(\Vec{r}) c(\Vec{r'}) \ \bigg| \ \frac{1}{\Vec{r}-\Vec{r'}} \ \bigg| \ b(\Vec{r}) d(\Vec{r'}) \ \bigg> \equiv \int&space;d&space;\Vec{r}&space;d&space;\Vec{r'}&space;\&space;\frac{&space;\&space;\psi^*_{n_a,l_a,m_a}&space;(\Vec{r})&space;\&space;\psi^*_{n_c,l_c,m_c}&space;(\Vec{r'})&space;\&space;\psi_{n_b,l_b,m_b}&space;(\Vec{r})&space;\&space;\psi_{n_d,l_d,m_d}&space;(\Vec{r'})&space;}{|\Vec{r}&space;-&space;\Vec{r'}|}" />

Where <img src="https://latex.codecogs.com/svg.latex?\psi_{n,l,m}(\Vec{r})" title="\psi_{n,l,m}(\Vec{r})" /> is a hydrogen-like atomic orbital with charge Z. These matrix elements occur, for example, in Hartree-Fock theory for a multi electron atom when a hydrogen-orbital basis is used.

The program is designed to be as simple as possible to use. For example, to generate all possible matrix elements up to principal quantum number n=2 and write them to a file called 'outfile.csv', type:

**python3 hydrogenic_matrix_elements.py -n_cutoff 2 -write_to_file outfile**

More examples can be found in the examples directory.

The only pre-requisite is a working installation of the [SymPy](https://www.sympy.org/en/index.html) library, which is used for the definitions of the <img src="https://latex.codecogs.com/svg.latex?\psi_{n,l,m}(\Vec{r})" title="\psi_{n,l,m}(\Vec{r})" /> and for its Wigner-3j symbol function.

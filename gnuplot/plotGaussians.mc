/* vim: ft=maxima
*/

load(numericalio);

width: 0.03;
min: 0.2;
max: 0.8;
num: 21;
kT: 2.494339;
filename: "ord20_s003_tdwt/5/coeffs.final";


gaussian(x,mu,sigma) := exp((-0.5)*((x-mu)/sigma)^2);
mat: read_matrix(filename);
positions: makelist(min+((i)*(max-min))/num,i,0,num-1);
allgauss(x) := args(rest(transpose(mat)[2],1))*gaussian(x,positions,width)/kT;
gausssum(x) := apply("+",allgauss(x));
colors: append([color],append([red],makelist(blue,i,0,num-1)));
plot2d(
  append([gausssum(x)],allgauss(x)),[x,min,max],
  [legend,false],
  colors,
  [xlabel, "Distance / nm"],
  [ylabel, "Bias / kT"],
  [y, -11, 0.8],
  [title, "Gaussian, N=20, σ=0.03"],
  /*[pdf_file, "~/tex/BF_Gaussian/Images/sum_ord20_s003_tdwt_new.pdf"]*/
  [pdf_file, "~/sum_gauss.pdf"]
);

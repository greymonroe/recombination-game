# Linkage & Recombination Lab

An interactive teaching game for PLS 152 on genetic linkage and recombination. It runs in any browser, on phones and laptops.

**Play it: https://greymonroe.github.io/recombination-game/**

## What students do

The organism has two chromosomes (n = 2) carrying seven marker sites: A–D on chromosome 1 and E–G on chromosome 2. Each site is a SNP where the blue parent and the red parent carry different DNA bases. The two inbred parents are crossed, and the F1 is selfed.

- **Simulate offspring** makes one F2 plant and animates it: meiosis in the F1 for the egg and the pollen with crossovers shown, one chromatid of each chromosome going into each gamete, the gametes coming together, and the offspring's base calls placed under each marker on its chromosomes.
- **Simulate [x] offspring** makes many at once. Each offspring is drawn as chromosomes colored by parental origin (or by zygosity), so crossovers are visible.
- A heat map and a plot show the recombination frequency between every pair of markers, estimated from the offspring genotypes. As the F2 population grows the estimates settle: close markers recombine rarely, distant markers more often, and markers on different chromosomes level off at 50%.

## How it works

- Simulated crossovers fall at random along each chromosome with no interference (Haldane model). Map positions are in `app.js` (`CHRS`).
- Recombination frequencies are maximum-likelihood estimates from the 3 × 3 table of F2 genotype counts for each marker pair (EM algorithm), capped at 50%.
- Plain HTML, CSS, and JavaScript with no build step and no dependencies. Open `index.html` to run it locally.

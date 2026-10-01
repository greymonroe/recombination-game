# Linkage & Recombination Lab

An interactive teaching game for PLS 152 on genetic linkage and recombination. It runs in any browser, on phones and laptops.

**Play it: https://greymonroe.github.io/recombination-game/**

## What students do

The organism has two chromosomes (n = 2) carrying seven SNP markers: A–D on chromosome 1 and E–G on chromosome 2. At each marker Parent 1 (blue) and Parent 2 (red) carry different DNA bases. The parents are crossed, and the F1 is testcrossed back to Parent 2, so each offspring's genotype shows which gamete the F1 made.

- **Simulate offspring** makes one testcross offspring and animates it: meiosis in the F1 with crossovers shown, one chromatid of each chromosome going into the gamete, the F1 gamete joining the tester gamete, and the offspring's base calls placed under each marker on its chromosomes.
- **Simulate [x] offspring** makes many at once. Each offspring is drawn as its two chromosomes colored by parental origin (or by zygosity), so crossovers are visible.
- A heat map and a plot show the recombination frequency between every pair of markers. As the population grows the estimates settle: close markers recombine rarely, distant markers more often, and markers on different chromosomes level off at 50%.

## How it works

- Simulated crossovers fall at random along each chromosome with no interference (Haldane model). Map positions and SNP bases are in `app.js` (`CHRS`).
- Recombination frequency for a marker pair is recombinant offspring ÷ total offspring, counted directly from the testcross genotypes.
- Plain HTML, CSS, and JavaScript with no build step and no dependencies. Open `index.html` to run it locally.

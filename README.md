# Linkage & Recombination Lab

An interactive teaching game for PLS 152 on genetic linkage and recombination. It runs in any browser, on phones and laptops.

**Play it: https://greymonroe.github.io/recombination-game/**

## What students do

The organism has two chromosomes (n = 2) carrying seven markers: A–D on chromosome 1 and E–G on chromosome 2. Two inbred parents (blue and orange) are crossed, and the F1 is selfed.

1. **Build offspring by hand.** Walk through meiosis in the F1: tap the paired chromosomes to place crossovers, divide the cell, and pick one of the four gametes as the egg. Repeat for the pollen. The offspring shows its genotype at every marker. After a few offspring, a heat map shows the recombination frequency between every pair of markers, estimated from the offspring genotypes.
2. **Simulate a population.** Simulate one offspring at a time or hundreds at once. Each offspring is drawn as chromosomes colored by parental origin (or by zygosity), so crossovers are visible. As the F2 population grows the estimates settle: close markers recombine rarely, distant markers more often, and markers on different chromosomes level off at 50%.

## How it works

- Simulated crossovers fall at random along each chromosome with no interference (Haldane model). Map positions are in `app.js` (`CHRS`).
- Recombination frequencies are maximum-likelihood estimates from the 3 × 3 table of F2 genotype counts for each marker pair (EM algorithm), capped at 50%.
- Plain HTML, CSS, and JavaScript with no build step and no dependencies. Open `index.html` to run it locally.

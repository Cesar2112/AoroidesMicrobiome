# The holobiont *Agelas oroides*




[![Sesam_img](https://www.imbe.fr/local/cache-vignettes/L400xH186/d46a112bebd61c35-0c5b6.png?1668533164)](https://sesam-anr.imbe.fr/)


**Sponges** (phylum Porifera) are functionally important members of marine benthic communities. Having arisen some 600 million years ago, sponges are among the most ancient of all extant animal lineages, evolving in a sea of mircoorganism, potential pathogens and parasites. Sponges are well known to produce bioactive compounds. Since the beginning of the 1970s, many research groups around the world have carried out chemical investigation on *Agelas* spp., resulting in fruitful achievements. Their studies revealed that *Agelas* sponges harbor many bioactive secondary metabolites, including alkaloids (especially bromopyrrole derivatives), terpenoids, glycosphingolipids, carotenoids, fatty acids and meroterpenoids. This bioactive molecules, can be produced by the sponge itself, by their associated microorganism or by both.
The concept of sponge holobiont gathers the host animal, its associated microorganisms and the biological processed occurring between this entities. In fact, microrganism can constitute up to 35% of the total sponge biomass. 
Our objective is to describe the *Agelas oroides* holobiont by describing their cytological components and the microorganism associated with this specie. In the present repository we present the bioinformatic pipeline used to describe the microbial community associated to the sponge *Agelas oroides*.

![Agelas_img](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/Aoro_git.JPG)


**Figure 1.** *A. oroides* is a Mediterranean sponge inhabiting coraligenous formations to mesophotic depths (10-120 m).


## Sampling and metadata associated to analysis
A total of 6 specimens of *A. oroides* were collected in the Marseille region in two semi-dark underwater habitats, at depth between 15-20 m: “Tiboulen du Frioul” and “Grotte à Pérès”. Using a knife, sponges were totally detached from the substrate. Specimens were transported in plastic containers filled with sea water and placed in aquarium with running sea water for at least 24 h acclimation before manipulation. Aquaria temperature ranged from 15 to 20 °C with 38 ppm salinity. For all specimens, 3 fragments of about 1 cm^3 were cut using a sterilized scalpel depending on the body part of the sponge: the first fragment identified as the external part, is the outer layer (about 1 mm depth) with a dark orange coloration. The second fragment identified as the internal part has a lighter orange color with a rough consistency and corresponds to ~ 98% of the total sponge volume. The last fragment is the entire sponge kept with both internal and external parts. 

A traceability information sheet: [Sample Sheet](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/SampleSheet_AO.csv) was also completed and contained for each sponge samples: the date of harvest, its geographical localization, a sample code for identification and discrimination of the body parts being analyzed (external part, internal part and entire sponge).

![AoroidesInternalSchema](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/FIG1_Git.jpg)


**Figure 2.** (A) *Agelas oroides in situ*. (B) Schematic representation of sampling: each specimen was subsampled depending on body parts: external part, internal part and the entire sponge. External and internal sections were used for cellular dissociation and enrichment (see Fig. 3), each fraction and/or entire sponge fragment was used for molecular and cytological analyses. Manipulation of samples were carried out at ~10 °C

In adittion we conduct a sponge cells / extracellular Prokaryote enrichement following two centrifugation forces as explained in the following schema:

![Fig3](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/Fig2_V11_Git.jpg)

**Figure 3.** Protocol for the process dissociation and enrichement of eukariotic sponge cells (C1) and extracellular prokaryotes (C2) fractions.


## Bioinformatic pipeline:

### VTAM python based data pre-treatment
[VTAM pre-treatment](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/Vtam_pipeline.md)

### Statistical analysis conducted in R language

[Statistical analysis](https://github.com/Cesar2112/AoroidesMicrobiome/blob/main/ScriptAoroidesMic.R)

**References:**
* Chu M-J, Li M, Ma H, Li P-L, Li G-Q (2022) Secondary metabolites from marine sponges of the genus Agelas : a comprehensive update insight on structural diversity and bioactivity. RSC Advances 12:7789–7820.
* Hentschel U, Piel J, Degnan SM, Taylor MW (2012) Genomic insights into the marine sponge microbiome. Nat Rev Microbiol 10:641–654.
* Zhang H, Dong M, Chen J, Wang H, Tenney K, Crews P (2017) Bioactive secondary metabolites from the marine sponge genus Agelas. Marine Drugs 15:351.
* González A, Dubut V, Corse E, Mekdad R, Dechatre T, Castet U, Hebert R, Meglécz E (2023) VTAM: A robust pipeline for validating metabarcoding data using controls. Computational and Structural Biotechnology Journal 21:1151–1156.





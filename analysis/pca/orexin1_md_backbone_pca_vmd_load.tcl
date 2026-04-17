mol new /home/jallow/compstrucbio/project-work/analysis/pca/orexin1_md_backbone_pca_average.gro type gro waitfor all
mol addfile /home/jallow/compstrucbio/project-work/analysis/pca/orexin1_md_backbone_pc1_filtered.xtc type xtc waitfor all
mol delrep 0 top
mol representation NewCartoon
mol color Name
mol selection all
mol material AOShiny
mol addrep top

mol new /home/jallow/compstrucbio/project-work/analysis/pca/orexin1_md_backbone_pca_average.gro type gro waitfor all
mol addfile /home/jallow/compstrucbio/project-work/analysis/pca/orexin1_md_backbone_pc1_pc2_filtered.xtc type xtc waitfor all
mol delrep 0 top
mol representation NewCartoon
mol color ColorID 1
mol selection all
mol material Transparent
mol addrep top

# Lmp_ML
Machin learning study of mass and heat transfer behavior of Graphene Aerogel (GA)

**Ownership:** Honglin Liu (MSE Ph.D. at Georgia Tech)
**Warning:** PLEASE DO NOT COPY AND USE WITHOUT AUTHOR'S PERMISSION !!!

## Application Background
- Membrane distillation
- Graphene aerogel
- Convective heat transfer
- Water molecules diffusivity

## Input
- Temperature  
- 3D matrix representation of graphene aerogel  
  (each element represents the number of atoms per unit space)

## Output
- Water vapor diffusion coefficient  
- Thermal conductivity  
- Morphological parameters:  
  - Density  
  - Porosity  
  - Tortuosity  

## Model
3D CNN + Transformer

## Input Data Source
Computed using LAMMPS:  
- `len` = Average sheet length of graphene  
- `sigma` = Virtual particle zero potential point  
- `temp` = Testing temperature  

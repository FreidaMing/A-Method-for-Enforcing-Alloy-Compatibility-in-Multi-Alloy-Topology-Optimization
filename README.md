# A Method for Enforcing Alloy Compatibility in Multi-Alloy Topology Optimization

## Overview

Metal additive manufacturing processes such as Directed Energy Deposition (DED), including Wire Arc Additive Manufacturing (WAAM), enable components with spatially varying material properties by combining multiple alloys in a single structure. This capability provides new opportunities to balance structural performance, material cost, and manufacturability.

However, conventional multi-material topology optimization (MMTO) methods often treat candidate materials as mechanically distinct but chemically independent. During fabrication, directly adjacent alloys may react and form brittle, weak, or otherwise undesirable phases. A mechanically efficient design may therefore be difficult or unsafe to manufacture when alloy compatibility is not considered.

This project develops a **Multi-Alloy Topology Optimization (MATO)** framework that explicitly incorporates alloy compatibility into the optimization process. The method combines:

- alloy-gradient design for identifying feasible composition pathways;
- graph-based pathfinding for selecting compatible alloy sequences;
- multi-material topology optimization for distributing alloys and material within the design domain; and
- compatibility and ordering constraints for preventing infeasible alloy adjacency.

The resulting framework simultaneously considers structural performance and metallurgical feasibility.

## Motivation

Existing MMTO methods can determine where different materials should be placed, but they often do not account for reactions between neighboring alloys during additive manufacturing. This omission can produce designs in which:

- chemically incompatible alloys are placed next to one another;
- harmful intermetallic phases may form during fabrication;
- abrupt composition changes reduce manufacturability; and
- the optimized design cannot be translated directly into a feasible deposition plan.

The proposed MATO framework addresses this gap by requiring material transitions to follow a predefined sequence of compatible alloys.

## Main Contributions

1. **Alloy-compatibility-aware topology optimization**  
   Alloy compatibility is incorporated directly into the optimization process rather than evaluated only after optimization.

2. **Integration of gradient-alloy design and structural optimization**  
   Feasible alloy transition pathways are identified and used to guide material placement.

3. **Enforcement of ordered alloy transitions**  
   Candidate alloys are arranged according to a compatible hierarchy so that nonadjacent or incompatible alloys are not placed in direct contact.

4. **Simultaneous structural and material design**  
   The method optimizes both structural topology and the spatial distribution of multiple alloys.

5. **Manufacturing-oriented design for DED and WAAM**  
   The optimized material layout is intended to remain consistent with the metallurgical requirements of multi-alloy additive manufacturing.

## Methodology

### 1. Compatible Alloy Path Identification

Candidate alloys are represented in a compatibility graph. Each node corresponds to an alloy, and each edge represents a feasible transition between two alloys.

A graph-based pathfinding procedure determines an ordered sequence of alloys connecting the desired material endpoints while avoiding incompatible transitions.

```text
Alloy 1 -> Alloy 2 -> Alloy 3 -> Alloy 4 -> Alloy 5
```

Only neighboring alloys in this sequence are allowed to form direct interfaces.

### 2. Multi-Alloy Topology Optimization

The structural design is represented using multiple level-set fields. The optimization minimizes a structural objective, such as compliance, while satisfying material-volume, ordering, and manufacturability constraints.


### 3. Alloy-Ordering Constraint

The material hierarchy is enforced so that higher-order alloys can only exist inside regions already occupied by the preceding alloy phase.

This nested representation prevents incompatible alloys from becoming direct neighbors and promotes gradual alloy transitions.

### 4. Reaction-Diffusion Level-Set Evolution

The level-set fields are updated using a reaction-diffusion equation.

## Computational Workflow

```text
Candidate alloy database
        |
        v
Compatibility assessment
        |
        v
Graph-based alloy pathfinding
        |
        v
Compatible alloy sequence
        |
        v
Initialization of level-set fields
        |
        v
Finite-element analysis
        |
        v
Sensitivity calculation
        |
        v
Reaction-diffusion update
        |
        v
Compatibility and hierarchy enforcement
        |
        v
Convergence check
        |
        v
Optimized topology and alloy distribution
```



## Running the Code

1. Clone the repository:

```bash
git clone https://github.com/FreidaMing/A-Method-for-Enforcing-Alloy-Compatibility-in-Multi-Alloy-Topology-Optimization.git
cd A-Method-for-Enforcing-Alloy-Compatibility-in-Multi-Alloy-Topology-Optimization
```

2. Open the project in MATLAB.

3. Review the model settings, including:

- design-domain dimensions;
- mesh resolution;
- boundary conditions;
- applied loads;
- material elastic properties;
- alloy cost or density;
- material-volume limits;
- compatible alloy order; and
- convergence parameters.

4. Run the main optimization script.

5. Review the generated outputs, including:

- structural topology;
- alloy distribution;
- level-set contours;
- objective and constraint histories;
- displacement field; and
- von Mises stress distribution.

## Expected Outputs

The framework can generate:

- optimized structural topology;
- spatial distribution of multiple alloys;
- compatible and ordered material transitions;
- compliance and volume histories;
- alloy-interface contours;
- displacement fields;
- stress distributions; and
- intermediate optimization results.

## Results

<img width="1203" alt="Optimized multi-alloy topology" src="https://github.com/user-attachments/assets/d1f56335-c4fb-426a-9b55-9d01cfbe943f">


## Applications

Potential applications include:

- functionally graded structural components;
- multi-alloy DED and WAAM;
- high-temperature structures;
- wear- or corrosion-resistant components;
- lightweight structures with locally tailored stiffness;
- thermal-management components; and
- mechanically optimized parts requiring metallurgically feasible transitions.

## Limitations and Ongoing Work

Current and future research directions include:

- strict control of minimum feature and layer thickness;
- improved manufacturing constraints;
- thermal and thermomechanical objectives;
- material cost as an objective or constraint;
- experimental validation of optimized alloy transitions;
- uncertainty quantification for alloy properties and compatibility data; and
- neural-network-based formulations for improved scalability.

## Publication

This repository supports the work:

> **A Method for Enforcing Alloy Compatibility in Multi-Alloy Topology Optimization**

The study was developed for the ASME International Design Engineering Technical Conferences and Computers and Information in Engineering Conference.

Replace the placeholders below with the final publication information:

```bibtex
@inproceedings{shu2025method,
  title={A Method for Enforcing Alloy Compatibility in Multi-Alloy Topology Optimization},
  author={Shu, Yalan and Allen, Marshall D and Vela, Brent and Malak, Richard},
  booktitle={International Design Engineering Technical Conferences and Computers and Information in Engineering Conference},
  volume={89237},
  pages={V03BT03A007},
  year={2025},
  organization={American Society of Mechanical Engineers}
}
```

## Contact

**Yalan Shu**  
Texas A&M University  
Email: `shuyl@tamu.edu`


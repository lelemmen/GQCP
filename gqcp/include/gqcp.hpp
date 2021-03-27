// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Basis/Integrals/BaseOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/Integrals/IntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/OneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/OneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/BaseMatrixPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/BaseScalarPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/BaseVectorPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/DoubleLondonHermiteCoulombIntegral.hpp"
#include "Basis/Integrals/Primitive/FunctionalOneElectronPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/FunctionalTwoElectronPrimitiveIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/HermiteCoulombIntegral.hpp"
#include "Basis/Integrals/Primitive/LondonHermiteCoulombIntegral.hpp"
#include "Basis/Integrals/Primitive/McMurchieDavidsonCoefficient.hpp"
#include "Basis/Integrals/Primitive/PrimitiveAngularMomentumIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveCanonicalKineticEnergyIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveCoulombRepulsionIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveElectronicDipoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveElectronicQuadrupoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveLinearMomentumIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveNuclearAttractionIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/Integrals/TwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/TwoElectronIntegralEngine.hpp"
#include "Basis/MullikenPartitioning/GMullikenPartitioning.hpp"
#include "Basis/MullikenPartitioning/RMullikenPartitioning.hpp"
#include "Basis/MullikenPartitioning/UMullikenPartitioning.hpp"
#include "Basis/MullikenPartitioning/UMullikenPartitioningComponent.hpp"
#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/OccupationType.hpp"
#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Basis/SpinorBasis/SpinResolvedOrbitalSpace.hpp"
#include "Basis/SpinorBasis/Spinor.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasisComponent.hpp"
#include "Basis/Transformations/BasisTransformable.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Basis/Transformations/JacobiRotatable.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "Basis/Transformations/OrbitalRotationGenerators.hpp"
#include "Basis/Transformations/RTransformation.hpp"
#include "Basis/Transformations/SimpleTransformation.hpp"
#include "Basis/Transformations/SpinResolvedBasisTransformable.hpp"
#include "Basis/Transformations/SpinResolvedJacobiRotatable.hpp"
#include "Basis/Transformations/UJacobiRotation.hpp"
#include "Basis/Transformations/UTransformation.hpp"
#include "Basis/Transformations/UTransformationComponent.hpp"
#include "Basis/Transformations/transform.hpp"
#include "DensityMatrix/DensityMatrixTraits.hpp"
#include "DensityMatrix/G1DM.hpp"
#include "DensityMatrix/G2DM.hpp"
#include "DensityMatrix/MixedSpinResolved2DMComponent.hpp"
#include "DensityMatrix/Orbital1DM.hpp"
#include "DensityMatrix/Orbital2DM.hpp"
#include "DensityMatrix/PureSpinResolved2DMComponent.hpp"
#include "DensityMatrix/Simple1DM.hpp"
#include "DensityMatrix/Simple2DM.hpp"
#include "DensityMatrix/SpinDensity1DM.hpp"
#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved1DMComponent.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"
#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/CompoundConvergenceCriterion.hpp"
#include "Mathematical/Algorithm/ConvergenceCriterion.hpp"
#include "Mathematical/Algorithm/FunctionalStep.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Algorithm/Step.hpp"
#include "Mathematical/Algorithm/StepCollection.hpp"
#include "Mathematical/Functions/BoysFunction.hpp"
#include "Mathematical/Functions/CartesianDirection.hpp"
#include "Mathematical/Functions/CartesianExponents.hpp"
#include "Mathematical/Functions/CartesianGTO.hpp"
#include "Mathematical/Functions/DyadicCartesianDirection.hpp"
#include "Mathematical/Functions/LinearCombination.hpp"
#include "Mathematical/Functions/ScalarFunction.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"
#include "Mathematical/Grid/CubicGrid.hpp"
#include "Mathematical/Grid/Field.hpp"
#include "Mathematical/Grid/WeightedGrid.hpp"
#include "Mathematical/Optimization/Accelerator/ConstantDamper.hpp"
#include "Mathematical/Optimization/Accelerator/DIIS.hpp"
#include "Mathematical/Optimization/ConsecutiveIteratesNormConvergence.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/CorrectionVectorCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/GuessVectorUpdate.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/MatrixVectorProductCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/ResidualVectorCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/ResidualVectorConvergence.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceMatrixCalculation.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceMatrixDiagonalization.hpp"
#include "Mathematical/Optimization/Eigenproblem/Davidson/SubspaceUpdate.hpp"
#include "Mathematical/Optimization/Eigenproblem/DenseDiagonalization.hpp"
#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "Mathematical/Optimization/LinearEquation/ColPivHouseholderQRSolution.hpp"
#include "Mathematical/Optimization/LinearEquation/HouseholderQRSolution.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"
#include "Mathematical/Optimization/Minimization/BaseHessianModifier.hpp"
#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "Mathematical/Optimization/Minimization/MinimizationEnvironment.hpp"
#include "Mathematical/Optimization/Minimization/Minimizer.hpp"
#include "Mathematical/Optimization/Minimization/NewtonStepUpdate.hpp"
#include "Mathematical/Optimization/Minimization/UnalteringHessianModifier.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NewtonStepUpdate.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"
#include "Mathematical/Optimization/NonLinearEquation/step.hpp"
#include "Mathematical/Optimization/OptimizationEnvironment.hpp"
#include "Mathematical/Representation/Array.hpp"
#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Mathematical/Representation/ImplicitMatrixSlice.hpp"
#include "Mathematical/Representation/ImplicitRankFourTensorSlice.hpp"
#include "Mathematical/Representation/LeviCivitaTensor.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Mathematical/Representation/MatrixRepresentationEvaluationContainer.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"
#include "Mathematical/Representation/StorageArray.hpp"
#include "Mathematical/Representation/Tensor.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Molecule/Nucleus.hpp"
#include "Molecule/elements.hpp"
#include "Operator/FirstQuantized/AngularMomentumOperator.hpp"
#include "Operator/FirstQuantized/BaseFQOperator.hpp"
#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"
#include "Operator/FirstQuantized/BaseReferenceDependentOperator.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/DiamagneticOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicQuadrupoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpin_zOperator.hpp"
#include "Operator/FirstQuantized/FQMolecularHamiltonian.hpp"
#include "Operator/FirstQuantized/FQMolecularMagneticHamiltonian.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/LinearMomentumOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"
#include "Operator/FirstQuantized/OrbitalZeemanOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/FirstQuantized/SpinZeemanOperator.hpp"
#include "Operator/SecondQuantized/EvaluatableScalarRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/AdjacencyMatrix.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"
#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"
#include "Operator/SecondQuantized/OperatorTraits.hpp"
#include "Operator/SecondQuantized/PureUSQTwoElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/SQOperatorStorage.hpp"
#include "Operator/SecondQuantized/SQOperatorStorageBase.hpp"
#include "Operator/SecondQuantized/SimpleSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SimpleSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "Physical/HomogeneousMagneticField.hpp"
#include "Physical/units.hpp"
#include "Processing/Properties/BaseElectricalResponseSolver.hpp"
#include "Processing/Properties/CIElectricalResponseSolver.hpp"
#include "Processing/Properties/DysonOrbital.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "Processing/Properties/properties.hpp"
#include "Processing/Properties/vAP1roGElectricalResponseSolver.hpp"
#include "QCMethod/CC/CCD.hpp"
#include "QCMethod/CC/CCDAmplitudesUpdate.hpp"
#include "QCMethod/CC/CCDEnergyCalculation.hpp"
#include "QCMethod/CC/CCDIntermediatesUpdate.hpp"
#include "QCMethod/CC/CCDSolver.hpp"
#include "QCMethod/CC/CCSD.hpp"
#include "QCMethod/CC/CCSDAmplitudesUpdate.hpp"
#include "QCMethod/CC/CCSDEnergyCalculation.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/CCSDIntermediatesUpdate.hpp"
#include "QCMethod/CC/CCSDSolver.hpp"
#include "QCMethod/CC/T2DIIS.hpp"
#include "QCMethod/CC/T2ErrorCalculation.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"
#include "QCMethod/Geminals/AP1roGJacobiOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/PSEnvironment.hpp"
#include "QCMethod/Geminals/vAP1roG.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFDensityMatrixCalculation.hpp"
#include "QCMethod/HF/GHF/GHFElectronicEnergyCalculation.hpp"
#include "QCMethod/HF/GHF/GHFErrorCalculation.hpp"
#include "QCMethod/HF/GHF/GHFFockMatrixCalculation.hpp"
#include "QCMethod/HF/GHF/GHFFockMatrixDIIS.hpp"
#include "QCMethod/HF/GHF/GHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/GHF/GHFSCFEnvironment.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFDensityMatrixCalculation.hpp"
#include "QCMethod/HF/RHF/RHFDensityMatrixDamper.hpp"
#include "QCMethod/HF/RHF/RHFElectronicEnergyCalculation.hpp"
#include "QCMethod/HF/RHF/RHFErrorCalculation.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixCalculation.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixDIIS.hpp"
#include "QCMethod/HF/RHF/RHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCMethod/HF/UHF/UHF.hpp"
#include "QCMethod/HF/UHF/UHFDensityMatrixCalculation.hpp"
#include "QCMethod/HF/UHF/UHFElectronicEnergyCalculation.hpp"
#include "QCMethod/HF/UHF/UHFErrorCalculation.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixCalculation.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixDIIS.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "QCMethod/HF/UHF/UHFSCFSolver.hpp"
#include "QCMethod/OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/JacobiOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"
#include "QCMethod/OrbitalOptimization/NewtonOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCMethod/RMP2/RMP2.hpp"
#include "QCModel/CC/CCD.hpp"
#include "QCModel/CC/CCSD.hpp"
#include "QCModel/CC/T1Amplitudes.hpp"
#include "QCModel/CC/T2Amplitudes.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "QCModel/Geminals/AP1roG.hpp"
#include "QCModel/Geminals/AP1roGGeminalCoefficients.hpp"
#include "QCModel/Geminals/APIGGeminalCoefficients.hpp"
#include "QCModel/Geminals/GeminalCoefficientsInterface.hpp"
#include "QCModel/HF/GHF.hpp"
#include "QCModel/HF/RHF.hpp"
#include "QCModel/HF/StabilityMatrices/GHFStabilityMatrices.hpp"
#include "QCModel/HF/StabilityMatrices/RHFStabilityMatrices.hpp"
#include "QCModel/HF/StabilityMatrices/UHFStabilityMatrices.hpp"
#include "QCModel/HF/UHF.hpp"
#include "QuantumChemical/DoublySpinResolvedBase.hpp"
#include "QuantumChemical/Spin.hpp"
#include "QuantumChemical/SpinResolved.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"
#include "QuantumChemical/spinor_tags.hpp"
#include "Utilities/CRTP.hpp"
#include "Utilities/Eigen.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/literals.hpp"
#include "Utilities/memory.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/type_traits.hpp"

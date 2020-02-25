//  MABE is a product of The Hintza Lab @ MSU
//     for general research information:
//         http://hintzelab.msu.edu/
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//          github.com/Hintzelab/MABE/wiki

//  This file was auto-generated with MBuilder.py

#ifndef __AutoBuild__Modules__
#define __AutoBuild__Modules__
#include "World/NKWorld/NKWorld.h"
#include "Genome/CircularGenome/CircularGenome.h"
#include "Brain/ConstantValuesBrain/ConstantValuesBrain.h"
#include "Optimizer/SimpleOptimizer/SimpleOptimizer.h"

#include "Archivist/DefaultArchivist.h"
#include "Archivist/LODwAPArchivist/LODwAPArchivist.h"
#include "Archivist/SSwDArchivist/SSwDArchivist.h"


//create a world
std::shared_ptr<AbstractWorld> makeWorld(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractWorld> newWorld;
  bool found = false;
  std::string worldType = AbstractWorld::worldTypePL->get(PT);
  if (worldType == "NK") {
    newWorld = std::make_shared<NKWorld>(PT);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find WORLD-worldType \"" << worldType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newWorld;
}


//create an optimizer
std::shared_ptr<AbstractOptimizer> makeOptimizer(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractOptimizer> newOptimizer;
  bool found = false;
  std::string optimizerType = AbstractOptimizer::Optimizer_MethodStrPL->get(PT);
  if (optimizerType == "Simple") {
    newOptimizer = std::make_shared<SimpleOptimizer>(PT);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find OPTIMIZER-optimizer \"" << optimizerType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newOptimizer;
}


//create an archivist
std::shared_ptr<DefaultArchivist> makeArchivist(std::vector<std::string> popFileColumns, std::shared_ptr<Abstract_MTree> _maxFormula, std::shared_ptr<ParametersTable> PT, std::string groupPrefix = ""){
  std::shared_ptr<DefaultArchivist> newArchivist;
  bool found = false;
  std::string archivistType = DefaultArchivist::Arch_outputMethodStrPL->get(PT);
  if (archivistType == "LODwAP") {
    newArchivist = std::make_shared<LODwAPArchivist>(popFileColumns, _maxFormula, PT, groupPrefix);
    found = true;
    }
  if (archivistType == "SSwD") {
    newArchivist = std::make_shared<SSwDArchivist>(popFileColumns, _maxFormula, PT, groupPrefix);
    found = true;
    }
  if (archivistType == "Default") {
    newArchivist = std::make_shared<DefaultArchivist>(popFileColumns, _maxFormula, PT, groupPrefix);
    found = true;
    }
  if (!found){
    std::cout << "  ERROR! could not find ARCHIVIST-outputMethod \"" << archivistType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newArchivist;
}


//create a template genome
std::shared_ptr<AbstractGenome> makeTemplateGenome(std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractGenome> newGenome;
  bool found = false;
  std::string genomeType = AbstractGenome::genomeTypeStrPL->get(PT);
  if (genomeType == "Circular") {
    newGenome = CircularGenome_genomeFactory(PT);
    found = true;
    }
  if (found == false){
    std::cout << "  ERROR! could not find GENOME-genomeType \"" << genomeType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newGenome;
}


//create a template brain
std::shared_ptr<AbstractBrain> makeTemplateBrain(int inputs, int outputs, std::shared_ptr<ParametersTable> PT){
  std::shared_ptr<AbstractBrain> newBrain;
  bool found = false;
  std::string brainType = AbstractBrain::brainTypeStrPL->get(PT);
  if (brainType == "ConstantValues") {
    newBrain = ConstantValuesBrain_brainFactory(inputs, outputs, PT);
    found = true;
    }
  if (found == false){
    std::cout << "  ERROR! could not find BRAIN-brainType \"" << brainType << "\".\n  Exiting." << std::endl;
    exit(1);
    }
  return newBrain;
}


//configure Defaults and Documentation
void configureDefaultsAndDocumentation(){
  Parameters::root->setParameter("BRAIN-brainType", (std::string)"ConstantValues");
  Parameters::root->setDocumentation("BRAIN-brainType", "brain to be used, [ConstantValues]");

  Parameters::root->setParameter("GENOME-genomeType", (std::string)"Circular");
  Parameters::root->setDocumentation("GENOME-genomeType", "genome to be used, [Circular]");

  Parameters::root->setParameter("ARCHIVIST-outputMethod", (std::string)"LODwAP");
  Parameters::root->setDocumentation("ARCHIVIST-outputMethod", "output method, [LODwAP, SSwD, Default]");

  Parameters::root->setParameter("OPTIMIZER-optimizer", (std::string)"Simple");
  Parameters::root->setDocumentation("OPTIMIZER-optimizer", "optimizer to be used, [Simple]");

  Parameters::root->setParameter("WORLD-worldType", (std::string)"NK");
  Parameters::root->setDocumentation("WORLD-worldType","world to be used, [NK]");
}


#endif /* __AutoBuild__Modules__ */

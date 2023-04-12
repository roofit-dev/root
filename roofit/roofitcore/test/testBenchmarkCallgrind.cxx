#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooFit/TestStatistics/buildLikelihood.h"
#include "RooAbsPdf.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "TFile.h"
#include "RooFit/MultiProcess/Config.h"
#include "RooFit/TestStatistics/RooRealL.h"
#include "RooMinimizer.h"

#include <memory>
#include <chrono>
#include <ctime>

#include <fstream>
#include <string>

#include "gtest/gtest.h"

TEST(TestBenchmarkCallgrind, testBenchmarkCallgrind)
{
	  auto start = std::chrono::system_clock::now();
    TFile *f = TFile::Open("/project/atlas/users/zwolffs/qt/roofit/benchmark_roofit/data/workspaces/2023/WS-Comb-STXS_toy.root");

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    RooWorkspace* w = (RooWorkspace*)f->Get("combWS");

    // Fixes for known features, binned likelihood optimization
    RooFIter iter = w->components().fwdIterator();
    RooAbsArg *arg;
    while ((arg = iter.next())) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
            std::cout << "Activating binned likelihood attribute for "
                                << arg->GetName() << std::endl;
        }
    }

    RooAbsData* data = w->data("toyData");
    auto mc = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
    auto global_observables = mc->GetGlobalObservables();
    auto nuisance_parameters = mc->GetNuisanceParameters();
    auto pdf = w->pdf(mc->GetPdf()->GetName());

    RooMinimizer* m;

    // Set queue type
    RooFit::MultiProcess::Config::Queue::setQueueType(RooFit::MultiProcess::Config::Queue::QueueType::Priority);

    // Parallel likelihood minimization
    // RooFit::MultiProcess::Config::setDefaultNWorkers(1);
    RooFit::MultiProcess::Config::LikelihoodJob::defaultNEventTasks = 1;
    RooFit::MultiProcess::Config::LikelihoodJob::defaultNComponentTasks = 1000;  // just a very large number, which I hope is larger than the actual number of components -> every component gets a task

    std::unique_ptr<RooAbsReal> likelihoodAbsReal{pdf->createNLL(
        *data, RooFit::Constrain(*nuisance_parameters),
        RooFit::GlobalObservables(*global_observables), RooFit::ModularL(true), RooFit::Offset(true))};

    RooMinimizer::Config cfg;
    cfg.parallelize = 1;
    cfg.timingAnalysis = false;
    cfg.enableParallelGradient = true;
    cfg.enableParallelDescent = true;
    cfg.offsetting = true;
    m = new RooMinimizer(*likelihoodAbsReal, cfg);
    printf("Modular style parallel likelihood minimization\n");
  
    m->setOffsetting(true);
    m->setMaxFunctionCalls(5650);
    m->setPrintLevel(1);
    m->setStrategy(0);
    m->setProfile(false);
    m->optimizeConst(2);
    m->setMinimizerType("Minuit2");
    m->setEps(1);

    // ------------------ MINIMIZATION ------------------ //

    m->minimize("Minuit2", "Migrad");

    delete m;
}
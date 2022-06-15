#include "RooFit.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
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

#include <gtest/gtest.h>

#include <RooTrace.h>

#include <fstream>
#include <string>

TEST(Bench, Bench)
{
    int nWorkers = 1;
    char const * fname = "/home/zef/Documents/phd/qt/roofit/benchmark_roofit/data/workspaces/HGam_mu.root";

    auto start = std::chrono::system_clock::now();
    TFile *f = TFile::Open(fname);

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

    RooAbsData* data = w->data("combData");
    auto mc = dynamic_cast<RooStats::ModelConfig *>(w->genobj("ModelConfig"));
    auto global_observables = mc->GetGlobalObservables();
    auto nuisance_parameters = mc->GetNuisanceParameters();
    auto pdf = w->pdf(mc->GetPdf()->GetName());

    // RooArgSet vars = w->allVars();

    // // Fixes for known features, binned likelihood optimization
    //  RooFIter iter_param = vars.fwdIterator();
    //  RooRealVar *arg_param;
    //  while ((arg_param = (RooRealVar*)iter_param.next())) {
    //      if (( ((TString)arg_param->GetName()).Contains("ATLAS") || ((TString)arg_param->GetName()).Contains("nbkg") ))
    //      {
    //          std::cout << "setting " << arg_param->GetName() << "const" << std::endl;
    //          arg_param->setConstant(1);
    //      }
    // }

    // ((RooRealVar*)w->arg("mu"))->setConstant(0);

    RooFit::MultiProcess::Config::setDefaultNWorkers(nWorkers);

    std::shared_ptr<RooFit::TestStatistics::RooAbsL> likelihood = RooFit::TestStatistics::buildLikelihood(pdf, data, RooFit::TestStatistics::ConstrainedParameters(*nuisance_parameters), RooFit::TestStatistics::GlobalObservables(*global_observables));
    RooMinimizer* m = new RooMinimizer(likelihood, RooFit::TestStatistics::LikelihoodMode::serial,
                   RooFit::TestStatistics::LikelihoodGradientMode::multiprocess);
    m->setOffsetting(true);

    m->setPrintLevel(1);
    m->setStrategy(0);
    m->setProfile(true);
    m->optimizeConst(2);
    m->setMinimizerType("Minuit2");
    //m->setVerbose(kTRUE);
    m->setEps(1);
    m->setEvalErrorWall(1); 

    long setup_t = std::chrono::duration_cast<std::chrono::milliseconds> (std::chrono::system_clock::now() - start).count();
    
    // ------------------ MINIMIZATION ------------------ //

    clock_t t_cpu = clock(); // CPU time
    start = std::chrono::system_clock::now(); // Wall time

    RooTrace::callgrind_zero();
    m->migrad();
    RooTrace::callgrind_dump();

    long minimize_t = std::chrono::duration_cast<std::chrono::milliseconds> (std::chrono::system_clock::now() - start).count(); // Wall time
    t_cpu = clock() - t_cpu; // CPU time

    m->save()->printValue(std::cout);

    m->cleanup(); // necessary in tests to clean up global _theFitter
    delete m;

    // ------------------ WRITE TIMING RESULTS TO FILE ------------------ //

    std::ofstream myfile;
    myfile.open("timings.csv", std::ios_base::app);
    myfile << std::to_string(setup_t) << "," << std::to_string(minimize_t) << "," << std::to_string(((float)t_cpu/CLOCKS_PER_SEC)*1000) << "\n";
    myfile.close();
}

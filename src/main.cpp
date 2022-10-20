#include <storm-config.h>
#include <boost/optional.hpp>

#include <storm-pars/api/storm-pars.h>
#include <storm-parsers/api/storm-parsers.h>
#include <storm/api/storm.h>

#include <storm-pars/api/region.h>

#include <storm-parsers/parser/PrismParser.h>

#include <storm/storage/prism/Program.h>
#include <storm/storage/jani/Property.h>

#include <storm/modelchecker/results/CheckResult.h>

#include <storm-cli-utilities/cli.h>
#include <storm-cli-utilities/model-handling.h>

#include <storm-pars/modelchecker/region/SparseParameterLiftingModelChecker.h>
#include <storm-pars/modelchecker/region/SparseDtmcParameterLiftingModelChecker.h>

#include <storm/models/ModelBase.h>

#include <storm/utility/initialize.h>
#include <storm/utility/Stopwatch.h>
#include <storm/utility/macros.h>
#include <storm/utility/Engine.h>
#include <storm-pars/utility/parameterlifting.h>
#include <storm-pars/utility/parametric.h>
#include "storm-pars/utility/ModelInstantiator.h"

#include <storm/adapters/RationalFunctionAdapter.h>

#include <storm-pars/transformer/SparseParametricDtmcSimplifier.h>

#include <storm/solver/stateelimination/NondeterministicModelStateEliminator.h>

#include "storm/models/sparse/StandardRewardModel.h"

//settings
#include <storm/settings/SettingsManager.h>
#include "storm-pars/settings/ParsSettings.h"
#include "storm-pars/settings/modules/ParametricSettings.h"
#include "storm-pars/settings/modules/RegionSettings.h"
//#include "storm-pars/settings/modules/MonotonicitySettings.h"

#include "storm/settings/modules/CoreSettings.h"
#include "storm/settings/modules/IOSettings.h"
#include "storm/settings/modules/BuildSettings.h"
#include "storm/settings/modules/ModelCheckerSettings.h"
#include "storm/settings/modules/EigenEquationSolverSettings.h"
#include "storm/settings/modules/EliminationSettings.h"

//inludes for other functions
#include <tuple>
#include <iostream>
#include <vector>
#include <cmath>
#include <boost/optional/optional_io.hpp>
#include <stdexcept>
#include <sstream>
#include <iomanip>

int numberOfRandom = 3;
boost::optional <uint_fast64_t> depthLimit = boost::none;
double refThreshold = 0.05;
bool faster = false;

/*!
 * Function to make string into double.
 *
 * @param str - String.
 *
 * @return The double of the string.
 */
double stringToDouble(std::string str){
    std::string::size_type positionOfNumerator = str.find("/");
    if(positionOfNumerator == std::string::npos){ //if it is not a fraction just convert and return
        return std::stod(str);
    }
    double lb = std::stod((str.substr(0,positionOfNumerator)));
    double ub = std::stod((str.substr(positionOfNumerator + 1)));

    double res = (lb/ub);

    const double multiplier = std::pow(10.0, 3);
    return std::ceil(res * multiplier) / multiplier; //return lb/ub; for not rounded
}

/*!
 * Helper - function to take the best value for a parameter.
 *
 * @param initial - The initial instantiation.
 * @param lower - The lower bound.
 * @param upper - The upper bound.
 *
 * @return The new parameter value.
 */
std::pair<std::string, double> calcNewBoundsForPar(std::pair<std::string, double> initial, double lower, double upper){
    double take;
    if((initial.second - upper) <= 0 && (initial.second - lower)>=0){
        take = initial.second; //take the initial value if it is in the interval
    }else if(abs(initial.second - upper) <= abs(initial.second - lower)){
        take = upper;
    }else{
        take = lower;
    }
    return std::make_pair(initial.first, take);
}

/*!
 * Helper - function to extract instantiation from region string.
 *
 * @param par - The string of a single parameter.
 * @param initial - The initial instantiation, not always needed.
 *
 * @return The instantiation.
 */
std::pair<std::string, double> getInstantiationFromParBounds(std::string par, std::vector<std::pair<std::string, double>> initial = {{"", 0}}){
    std::string res = "";
    //!Get the position of the lower bound
    std::string::size_type positionOfFirstRelation = par.find("<=");
    STORM_LOG_THROW(positionOfFirstRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the first number");

    //!Get the position of the upper bound
    std::string::size_type positionOfSecondRelation = par.find("<=", positionOfFirstRelation+2);
    STORM_LOG_THROW(positionOfSecondRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the parameter");

    //!Get the parameter as string and the bounds as doubles
    std::string parameter = par.substr(positionOfFirstRelation+2,positionOfSecondRelation-(positionOfFirstRelation+2));
    parameter.erase(std::remove_if (parameter.begin(), parameter.end(), ::isspace), parameter.end());

    double lb = stringToDouble((par.substr(0,positionOfFirstRelation)));
    double ub = stringToDouble((par.substr(positionOfSecondRelation+2)));

    std::vector<std::pair<std::string, double>> temp = {{"", 0}}; //todo do it with boolean, probably more simple
    if(lb-ub != 0 && initial !=  temp){
        std::pair<std::string, double> found = {"", 0};
        std::pair<std::string, double> check = {"", 0};
        for(auto const& inst: initial){
            if(inst.first == parameter){
                found = inst;
                break;
            }
        }
        if(check == found){
            throw std::invalid_argument("Error in getInstantiationFromParBounds");
        }
        return calcNewBoundsForPar(found, lb, ub);
    }

    return std::make_pair(parameter, ub);
}

/*!
 * Function to extract instantiation from region string.
 *
 * @param parameters - A string of the current parameter bounds.
 * @param region_string - The constant, that is used to change the bounds, calculated in calc_bound.
 *
 * @return The instantiation.
 */
std::vector<std::pair<std::string, double>> getInstantiationFromRegion(std::set<storm::RationalFunctionVariable> parameters , std::string region_string,  std::vector<std::pair<std::string, double>> initial = {{"", 0}}){
    std::vector<std::pair<std::string, double>> res;
    std::vector<std::string> parameterBoundaries;
    boost::split(parameterBoundaries, region_string, boost::is_any_of(","));
    for (auto const& parameter : parameterBoundaries){
        if (!std::all_of(parameter.begin(),parameter.end(), ::isspace)){ //skip this string if it only consists of space
            //!Appends to res the bounds of the parameter
            res.push_back(getInstantiationFromParBounds(parameter, initial));
        }
    }

    return res;
}

/*!
 * Function to calculate the bounds of a single parameter, to later append to the region string.
 *
 * @param par - A string of the current parameter bounds.
 * @param bound - The constant, that is used to change the bounds, calculated in calc_bound.
 *
 * @return String representing the new parameter bounds.
 */
std::string change_region(std::string par, double bound){
    std::string res = "";
    //!Get the position of the lower bound
    std::string::size_type positionOfFirstRelation = par.find("<=");
    STORM_LOG_THROW(positionOfFirstRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the first number");

    //!Get the position of the upper bound
    std::string::size_type positionOfSecondRelation = par.find("<=", positionOfFirstRelation+2);
    STORM_LOG_THROW(positionOfSecondRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the parameter");

    //!Get the parameter as string and the bounds as doubles
    std::string parameter = par.substr(positionOfFirstRelation+2,positionOfSecondRelation-(positionOfFirstRelation+2));
    parameter.erase(std::remove_if (parameter.begin(), parameter.end(), ::isspace), parameter.end());

    double lb = std::stod((par.substr(0,positionOfFirstRelation)));
    double ub = std::stod((par.substr(positionOfSecondRelation+2)));

    //!Calculate the new bounds
    lb = lb - bound;

    ub = ub + bound;

    //!Removes the trailing zeros of a double
    std::ostringstream oss;
    oss << std::setprecision(8) << std::noshowpoint << lb;
    std::string str_lb = oss.str();
    oss.str("");
    oss << std::setprecision(8) << std::noshowpoint << ub;
    std::string str_ub = oss.str();

    //!Build the string as needed for a region, e.g. "0.1 <= p <= 0.9"
    res.append(str_lb);
    res.append("<=");
    res.append(parameter);
    res.append("<=");
    res.append(str_ub);
    res.append(",");

    return res;
}

/*!
 * Builds the region to give as an argument for the PLA, using change_region.
 *
 * @param parameters - The parameters of the model.
 * @param region_string - The initial values of the parameters, in most cases lb=ub, e.g. "0.5 <= p <= 0.5".
 * @param epsilon - The constant, that is used to change the bounds, calculated in calc_bound.
 *
 * @return The region string to use as an argument for the PLA algorithm.
 */
std::string build_region(std::set<storm::RationalFunctionVariable> parameters , std::string region_string, double epsilon){
    std::string res = "";
    std::vector<std::string> parameterBoundaries;
    boost::split(parameterBoundaries, region_string, boost::is_any_of(","));
    for (auto const& parameter : parameterBoundaries){
        if (!std::all_of(parameter.begin(),parameter.end(), ::isspace)){ //skip this string if it only consists of space
            //!Appends to res the bounds of the parameter
            res.append(change_region(parameter, epsilon));
        }
    }
    //!Removes the comma of the last parameter bounds
    res.pop_back();

    std::cout << "Doing PLA on: "<< res <<std::endl;
    return res;
}

std::string change_parSpace(std::string par){
    std::string res = "";
    //!Get the position of the lower bound
    std::string::size_type positionOfFirstRelation = par.find("<=");
    STORM_LOG_THROW(positionOfFirstRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the first number");

    //!Get the position of the upper bound
    std::string::size_type positionOfSecondRelation = par.find("<=", positionOfFirstRelation+2);
    STORM_LOG_THROW(positionOfSecondRelation!=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the region" << par << " I could not find a '<=' after the parameter");

    //!Get the parameter as string and the bounds as doubles
    std::string parameter = par.substr(positionOfFirstRelation+2,positionOfSecondRelation-(positionOfFirstRelation+2));
    parameter.erase(std::remove_if (parameter.begin(), parameter.end(), ::isspace), parameter.end());

    //!Removes the trailing zeros of a double
    std::ostringstream oss;
    oss << std::setprecision(8) << std::noshowpoint << 0;
    std::string str_lb = oss.str();
    oss.str("");
    oss << std::setprecision(8) << std::noshowpoint << 1;
    std::string str_ub = oss.str();

    //!Build the string as needed for a region, e.g. "0.1 <= p <= 0.9"
    res.append(str_lb);
    res.append("<=");
    res.append(parameter);
    res.append("<=");
    res.append(str_ub);
    res.append(",");

    return res;
}

std::string build_parSpace (std::set<storm::RationalFunctionVariable> parameters , std::string region_string){
    std::string res = "";
    std::vector<std::string> parameterBoundaries;
    boost::split(parameterBoundaries, region_string, boost::is_any_of(","));
    for (auto const& parameter : parameterBoundaries){
        if (!std::all_of(parameter.begin(),parameter.end(), ::isspace)){ //skip this string if it only consists of space
            //!Appends to res the bounds of the parameter
            res.append(change_parSpace(parameter));
        }
    }
    //!Removes the comma of the last parameter bounds
    res.pop_back();

    std::cout << "Doing PLA on: "<< res <<std::endl;
    return res;
}

/*!
 * Calculates the distance between two parametric instantiations.
 *
 * @param inst_old - The old instantiation of the parameters.
 * @param inst_new - The new instantiation of the parameters.
 * @param param_number - The number of parameters in the model.
 *
 * @return The EC-distance, with two decimal numbers.
 */
double EC_dist(std::vector<std::pair<std::string, double>> oldInst, std::vector<std::pair<std::string, double>> newInst){
    double temp = 0;
    for (int i = 0; i < newInst.size(); i++){
        for(int j = 0; j < newInst.size(); j++){
            if(newInst[i].first == oldInst[j].first){
                temp +=  std::pow((newInst[i].second - oldInst[j].second), 2); //(i_new - i_old)*(i_new - i_old)
                break;
            }
        }
    }
    double distance = sqrt(temp);

    //!Rounds the distance to 2 decimal places
    const double multiplier = std::pow(10.0, 3);
    return std::ceil(distance * multiplier) / multiplier; //return sqrt(temp); for not rounded
}

/*!
 * Calculates the constant to be used to change the parameter bounds in change_region.
 *
 * @param epsilon - The initial value given by the user.
 * @param numberOfParameters - The number of parameters in the model.
 *
 * @return The rounded to two decimal places constant.
 */
double calc_bound(double epsilon, size_t numberOfParameters){
    //!Checks if numberOfParameters can be safely casted to integer.
    if ( numberOfParameters > INT_MAX ){
        throw std::overflow_error("The number of parameters is to big" + numberOfParameters);
    }
    int parNumber = static_cast<int>(numberOfParameters);
    double res = epsilon/(sqrt(parNumber));

    //!Rounds to three decimal places
    const double multiplier = std::pow(10.0, 3);
    return std::ceil(res * multiplier) / multiplier; //return res; for not rounded
}

/*!
 * Gets a random region from the result vector of PLA, after only the satisfying regions are left.
 *
 * @param result_onlySAT - the vector with the satisfying regions.
 *
 * @return a random region-string.
 */
std::string getRandomRegion(std::vector<std::pair<storm::storage::ParameterRegion<storm::RationalFunction>, storm::modelchecker::RegionResult>> result_onlySAT){
    std::string region_string;
    int random = rand() % result_onlySAT.size();
    //! Get a random region
    auto const& region = result_onlySAT[random];

    //! Extract only the string of this region.
    std::stringstream first;
    first << region.first;
    region_string = first.str();
    return region_string;
}

/*!
 * Function of processing the CenterSat regions. It gets the center points of the regions and searches for the minimal one.
 *
 * @param result_exists - the vector with the CenterSat regions.
 * @param modelParameters - the parameters.
 * @param region_string - the initial instantiation, used for computing the distance.
 *
 * @return int - success 1 or failure 0.
 */
int processCenterSat(std::vector<std::pair<storm::storage::ParameterRegion<storm::RationalFunction>, storm::modelchecker::RegionResult>> result_exists, std::set<storm::RationalFunctionVariable> modelParameters, std::string region_string){
    if(result_exists.size()< numberOfRandom){
        std::cout << "There are only this much CenterSat regions:"<< result_exists.size() <<std::endl;
        numberOfRandom = result_exists.size();
    }
    std::vector<std::pair<std::string, double>> initialInst = getInstantiationFromRegion(modelParameters, region_string);
    std::vector<std::pair<std::string, double>> newInst;
    double distanceMin = 5;
    for (auto const& res : result_exists) {
        //todo remove later
        std::cout << res.second << ": " << res.first <<std::endl;

        std::ostringstream tempor;
        tempor << res.first.getCenterPoint();
        std::string val = tempor.str();
        val.erase(std::remove(val.begin(), val.end(), '{'), val.end());
        val.erase(std::remove(val.begin(), val.end(), '}'), val.end());

        std::vector<std::string> parameters;
        boost::split(parameters, val, boost::is_any_of(","));
        std::vector<std::pair<std::string, double>> candidate;

        for(auto const& par : parameters ){
            std::string::size_type positionOfFirstRelation = par.find(":");
            STORM_LOG_THROW(positionOfFirstRelation !=std::string::npos, storm::exceptions::InvalidArgumentException, "When parsing the valuation" << par << " I could not find a '='.");

            std::string parameter = par.substr(0, positionOfFirstRelation);
            parameter.erase(std::remove_if(parameter.begin(), parameter.end(), ::isspace), parameter.end());

            double value = stringToDouble((par.substr(positionOfFirstRelation + 1)));

            candidate.push_back(std::make_pair(parameter, value));
        }

        double candDistance = EC_dist(initialInst, candidate);
        if(candDistance < distanceMin){
            newInst = candidate;
            distanceMin = candDistance;
        }
        numberOfRandom--;
        if(numberOfRandom == 0){
            break;
        }
    }

    std::cout << distanceMin << " That is the newfound distance ;D" <<std::endl; //TODO fix that later

    //!Prints the instantation
    for (auto const& temp : newInst) {
        std::cout << temp.first << "= " << temp.second <<std::endl;
    }

    return 1;
}

/*!
 * Calculates the constant to be used to change the parameter bounds in change_region.
 *
 * @param path_to_model - The path to the model file, given by the user.
 * @param property_string - The property for which PLA would be executed.
 * @param region_string - The initial values of the parameters, in most cases lb=ub, e.g. "0.5 <= p <= 0.5".
 * @param epsilon - The initial value given by the user.
 * @param depthLimit - How many times a region would be split by PLA before declared for Unknown.
 * @param refThreshold - How big percent of the region can be left undecided by PLA.
 *
 * @return A vector containing only the SAT regions after the use of PLA.
 */
int pla(std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> model, std::vector<std::shared_ptr<const storm::logic::Formula>> formulae, std::string region_string, double epsilon, bool flag = false) {
     //auto model = storm::parser::DirectEncodingParser<storm::RationalFunction>::parseModel(path_to_model)->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>();

     //!Extract parameters
     auto modelParameters = storm::models::sparse::getProbabilityParameters(*model);
     auto rewParameters = storm::models::sparse::getRewardParameters(*model);
     modelParameters.insert(rewParameters.begin(), rewParameters.end());

     std::cout << modelParameters << std::endl;

     //! Uses the parameters and epsilon to calculate the constant to change the region and builds the new region
     std::string input;
     if (flag) {
        input = build_parSpace(modelParameters, region_string);
     }else{
         double bound = calc_bound(epsilon, modelParameters.size());
         input = build_region(modelParameters, region_string, bound);
     }
     auto region = storm::api::parseRegions<storm::RationalFunction>(input, *model);
    //! Preparations for Model Checker
    auto regionSettings = storm::settings::getModule<storm::settings::modules::RegionSettings>();
    auto engine = regionSettings.getRegionCheckEngine();
    boost::optional <uint_fast64_t> optionalDepthLimit = depthLimit;
    auto refinementThreshold = storm::utility::convertNumber<storm::RationalFunction>(refThreshold);
    //! Not using monotonicity. Comment this if you wan to use later version of storm.
    uint_fast64_t monThresh = 0;
    storm::api::MonotonicitySetting monotonicitySettings = storm::api::MonotonicitySetting();

    storm::utility::Stopwatch PLAwatch(true);
    try {
        auto result = storm::api::checkAndRefineRegionWithSparseEngine<storm::RationalFunction>(
                model, storm::api::createTask<storm::RationalFunction>(formulae[0], true), region.front(), engine,
                refinementThreshold,
                optionalDepthLimit, regionSettings.getHypothesis(), false);

    PLAwatch.stop();
    STORM_PRINT("Time for PLA: " << PLAwatch << ".\n\n");

    //! Extract only the SAT regions
    std::vector<std::pair<storm::storage::ParameterRegion<storm::RationalFunction>, storm::modelchecker::RegionResult>> result_onlySAT;
     std::vector<std::pair<storm::storage::ParameterRegion<storm::RationalFunction>, storm::modelchecker::RegionResult>> result_exists;
    for (auto const& res : result->getRegionResults()) {
        if (res.second == storm::modelchecker::RegionResult::AllSat) {
            result_onlySAT.push_back(res);
        }
    }

    if(result_onlySAT.empty()) {
        for (auto const &res: result->getRegionResults()) {
            if (res.second == storm::modelchecker::RegionResult::CenterSat) {
                result_exists.push_back(res);
            }
        }
		if(result_exists.empty()){
            if(flag){
                return 23;
            }
			return 9;
		}
        return processCenterSat(result_exists, modelParameters, region_string);
    }

    if(result_onlySAT.empty() && result_exists.empty()){
        std::cout << "No feasible solution found for this arguments. Try again." <<std::endl;
        return 9;
    }

    if(result_onlySAT.size()< numberOfRandom){
        std::cout << "There are only this much AllSatregions:"<< result_onlySAT.size() <<std::endl;
        numberOfRandom = result_onlySAT.size();
    }

    std::string reg = getRandomRegion(result_onlySAT);
    //result_onlySAT.remove_if(result_onlySAT.begin(),result_onlySAT.end(), isUsed(result_onlySAT));

    std::vector<std::pair<std::string, double>> initialInst = getInstantiationFromRegion(modelParameters, region_string);

    std::vector<std::pair<std::string, double>> newInst = getInstantiationFromRegion(modelParameters, reg, initialInst);

    double distanceMin = EC_dist(initialInst, newInst);

    while((numberOfRandom-1) > 0){
        if(result_onlySAT.empty()){
            //Should not happen because of the above check, but put it just to be sure
            std::cout << "The result vector is empty, yet you still try to get elements from it!" <<std::endl;
            break;
        }
        std::string newReg = getRandomRegion(result_onlySAT);
        std::vector<std::pair<std::string, double>> candidate = getInstantiationFromRegion(modelParameters, newReg, initialInst);
        double candDistance = EC_dist(initialInst, candidate);
        if(candDistance < distanceMin){
            newInst = candidate;
            distanceMin = candDistance;
        }
        numberOfRandom--;
    }


    std::cout << distanceMin << " The EC - distance between the original and new instantiation:" <<std::endl; 

    //!Prints the instantation
    for (auto const& temp : newInst) {
        std::cout << temp.first << "= " << temp.second <<std::endl;
    }

    }catch(storm::exceptions::InvalidArgumentException e){
        throw std::invalid_argument("Epsilon is to big and the bounded region is making the transition matrix non-stochastic. Please choose smaller epsilon or change initial values.");
    }

    return 1;
}

/*!
 * Sets the options from the user input.
 *
 * @param argc - The number of arguments given by the user.
 * @param argv - The arguments given by the user.
 */
void setOptions(int argc ,char *argv[]){
     for(int i = 5; i < argc; i++) {
         if (std::string(argv[i]) == "--numberOfRandom") {
             numberOfRandom = std::stoi(std::string(argv[i+1]) );
             std::cout << "Number of random instantiations: " << argv[i + 1] << endl;
         } else if (std::string(argv[i])  == "--depthLimit") {
             depthLimit = std::stoi(std::string(argv[i+1]));
             std::cout << "Depth limit: " << argv[i + 1] << endl;
         } else if (std::string(argv[i])  == "--refThreshold") {
             refThreshold = std::stod(std::string(argv[i+1]));
             if(refThreshold <= 0 || refThreshold >= 1){
                 std::cout << "Please, choose a suitable value for refThreshold in (0,1)." << endl;
             }
             std::cout << "Refinement threshold: " << argv[i + 1] << endl;
         } else if (std::string(argv[i])  == "--faster") {
             faster = true;
             std::cout << "Faster computation: " << argv[i + 1] << endl;
         }else{
             if(i%2 == 1){
                 std::cout << "There is no such option: " << argv[i] << endl;
             }
         }
     }
 }

/*!
 * Checks if the arguments are correct to start the algorithm.
 *
 * @param argc - The number of arguments given by the user.
 * @param argv - The arguments given by the user.
 */
void check_argv(int argc ,char *argv[]){
    std::ifstream test(argv[1]);
    if (!test){
        std::cout << argv[1]  << "is not a valid path/file." << endl;
        throw std::invalid_argument("");
    }
    std::cout << "Region: " << argv[3]  << endl;
    std::cout << "Epsilon: " << argv[4]  << endl;
    if(argc > 5){
        setOptions(argc, argv);
    }
}

/*!
 * argv[1] - file of the model
 * argv[2] - file of PCTL formula
 * argv[3] - instantiation in form of a region, e.g. "0.5 <= p <=0.5, 0.5 <= q <=0.5"
 * argv[4] - epsilon
 */
int main(int argc, char *argv[]) {
    //! Check number of arguments
    if(argc < 5){
        throw std::invalid_argument("The following arguments are obligatory: model, pctl, region's string and epsilon!");
    }
    //! Check if the arguments are correct
    check_argv(argc,argv);

    //! Init loggers.
    storm::utility::setUp();
    //! Init settings.
    storm::settings::initializeParsSettings("Storm-pars", "storm-pars");

    storm::utility::Stopwatch totalTimer(true);

    storm::utility::Stopwatch modelParsingWatch(true);

    std::string path_to_model = argv[1];
    std::string property_string = argv[2];

    //drn Parser
    auto model = storm::api::buildExplicitDRNModel<storm::RationalFunction>(argv[1])->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
    std::string formulasString = argv[2];
    auto formulae = storm::api::extractFormulasFromProperties(storm::api::parseProperties(formulasString));

    //prism Parser
    /*
    auto program = storm::parser::PrismParser::parse(path_to_model);
    assert(program.getModelType() == storm::prism::Program::ModelType::DTMC);
    auto properties = storm::api::parsePropertiesForPrismProgram(property_string, program);
    std::vector<std::shared_ptr<const storm::logic::Formula>> formulae = storm::api::extractFormulasFromProperties(properties);
    std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> model = storm::api::buildSparseModel<storm::RationalFunction>(program, formulae)->template as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
     */
    modelParsingWatch.stop();
    STORM_PRINT("Time for model input parsing: " << modelParsingWatch << ".\n\n");

    int res = pla(model, formulae, argv[3], std::stod(argv[4]));
    int counter = 1;

    //! Repeat till a feasible solution is found for epsilon, up to 5 times(ensures termination)
    while(res == 9 && counter <= 5){
        if(counter == 2 && faster == true) {
            int temp = pla(model, formulae, argv[3], std::stod(argv[4]) + 0.1 * counter, true);
            if (temp == 23) {
                std::cout << "No feasible solution exists for this property." << endl;
                return 0;
            }
            counter++;
        }else {
            res = pla(model, formulae, argv[3], std::stod(argv[4]) + 0.1 * counter);
            counter++;
        }
    }
    if(res == 9){
        std::cout << "No feasible solution found after increasing epsilon with 0.5 please give different arguments and try again." << argv[3]  << endl;
    }

    totalTimer.stop();
    storm::cli::printTimeAndMemoryStatistics(totalTimer.getTimeInMilliseconds());

    //!Clean after ready
    storm::utility::cleanUp();
    return 0;
}
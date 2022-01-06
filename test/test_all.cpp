/*
int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/

#include "../include/kde1d.hpp"

int main() {
    
    Eigen::VectorXd x(100);
    kde1d::Kde1d fit;
    fit.fit(x);
    fit.pdf(x);
    fit.cdf(x);
    fit.quantile(x.cwiseMax(0));

    std::cout << "success" << std::endl;

    return 0;
}

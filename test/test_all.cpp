/*
int main(int argc, char **argv) {

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/

#include "../include/kde1d.hpp"

int
main()
{

    {
        // continuous data
        Eigen::VectorXd x(100);
        kde1d::Kde1d fit;
        fit.fit(x);
        fit.pdf(x);
        fit.cdf(x);
        fit.quantile(x.cwiseMax(0));
    }

    {
        // discrete data
        Eigen::VectorXi x = Eigen::VectorXi::LinSpaced(101, -50, 50);
        kde1d::Kde1d fit;
        fit.fit(x);
        fit.pdf(x);
        fit.cdf(x);
        fit.quantile(Eigen::VectorXd::LinSpaced(100, 0.001, 0.999));
    }

    std::cout << "success" << std::endl;

    return 0;
}

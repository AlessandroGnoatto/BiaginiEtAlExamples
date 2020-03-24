package discva;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariableFromDoubleArray;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloBlackScholesModel;
import net.finmath.plots.Named;
import net.finmath.plots.Plot2D;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

public class DiscVAExample {

	public static void main(String[] args) throws CalculationException, IOException {
		double deltaT = 0.001;
		int numberOfTimeSteps = 1000;
		
		/*
		 * Creiamo una partizione del tempo su cui il processo verrà discretizzato.
		 */
		TimeDiscretization timeDiscretization =
				new TimeDiscretizationFromArray(0.0, numberOfTimeSteps, deltaT);
		
		/*
		 * Creiamo un modello finanziario (Il modello lognormale visto a lezione)
		 */
		int numberOfPaths = 10000;
		double initialValue = 100;
		double collateralRate = 0.01;
		double volatility = 0.25;
		
		MonteCarloBlackScholesModel bsModel = new MonteCarloBlackScholesModel(
				timeDiscretization,
				numberOfPaths,
				initialValue,
				collateralRate,
				volatility);
		
		
		//We price a forward.
		double strike =  80;
		double maturity = deltaT *(numberOfTimeSteps -1);
		int quantity = 1000;
		
		//V = S0 - K\exp(-rT)
		double exactPrice = quantity*(initialValue - strike * Math.exp(-collateralRate * maturity));
		double mcPrice = quantity* (Math.exp(-collateralRate * maturity)*(bsModel.getAssetValue(maturity, 0).getAverage() - strike));
		
		//PV from the trading desk, the deal is uncollateralized
		double unsecuredRate = 0.05;
		double frontOfficePrice = quantity * (initialValue*Math.exp(-(unsecuredRate - collateralRate) * maturity) - strike * Math.exp(-unsecuredRate * maturity));
		double frontOfficeMcPrice = quantity * (Math.exp(-unsecuredRate * maturity)*(bsModel.getAssetValue(maturity, 0).getAverage() - strike) );
		
		/*
		 * Compute the pathwise exposure according to the xVA desk
		 */
		RandomVariable[] pathwiseExposure = new RandomVariable[numberOfTimeSteps];
		pathwiseExposure[0] = new RandomVariableFromDoubleArray(timeDiscretization.getTime(0), exactPrice);
		
		for(int i = 0; i < numberOfTimeSteps; i++) {
			double runningTime = timeDiscretization.getTime(i);
			
			pathwiseExposure[i] = (bsModel.getAssetValue(timeDiscretization.getTime(i), 0)
					.mult(Math.exp(-collateralRate * (maturity - runningTime)))
					.sub(strike * Math.exp(-collateralRate * (maturity - runningTime)))).mult(quantity);
					
		}
		
		double discVA = 0.0;
		double spread = collateralRate - unsecuredRate;
		/*
		 * Compute the DiscVA.
		 */
		for(int i = 0; i < numberOfTimeSteps; i++) {
			double runningTime = timeDiscretization.getTime(i);
			//double runningTimePlusOne = timeDiscretization.getTime(i+1);
			discVA += pathwiseExposure[i].getAverage() * Math.exp(-unsecuredRate * runningTime);
						
		}
		
		discVA = discVA * spread * deltaT;
		
		System.out.println("xVA Desk analytical price: " + exactPrice);
		System.out.println("xVA Desk MC price: " + mcPrice);
		System.out.println("Front Office price: " + frontOfficePrice);
		System.out.println("Front Office MC price: " + frontOfficeMcPrice);
		System.out.println("Delta Mc Prices: "+ (frontOfficeMcPrice - mcPrice));
		System.out.println("The DiscVA is:");
		System.out.println(discVA);
		//System.out.println("Delta NPV between the desks" + (exactPrice - frontOfficePrice));
		System.out.println("Reconstruct xVA desk price from front Office price and DiscVA");
		System.out.println((frontOfficePrice - discVA));
		
		/*
		 * Just some nice plots...
		 */
		
		 DoubleUnaryOperator EPE = (time) -> {
		      int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
			  RandomVariable valuesPositiveExposure = pathwiseExposure[index].floor(0.0);

			  double expectedPositiveExposure   = valuesPositiveExposure.getAverage();

			  return expectedPositiveExposure;


		    };
		    
		    DoubleUnaryOperator ENE = (time) -> {
		    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
		    	RandomVariable valuesNegativeExposure = pathwiseExposure[index].cap(0.0);

		    	double expectedNegativeExposure   = valuesNegativeExposure.getAverage();

		    	return expectedNegativeExposure;

		    };
		    
		    DoubleUnaryOperator PFE = (time) -> {
		    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
		    	RandomVariable valuesEstimatedExposure;
		    	valuesEstimatedExposure = pathwiseExposure[index];

		    	double exposureQuantile = valuesEstimatedExposure.getQuantile(0.95);

		    	return exposureQuantile;


		    };
		    
		    Plot2D plot = new Plot2D(0.0, 0.99, 900,Arrays.asList(
		        new Named<DoubleUnaryOperator>("EPE", EPE),
		        new Named<DoubleUnaryOperator>("ENE", ENE),
		        new Named<DoubleUnaryOperator>("95%PFE", PFE)));
		    plot.setTitle("Expected Positive Exposure").setXAxisLabel("time").setYAxisLabel("Exposure").setIsLegendVisible(true);
		    plot.show();
		    plot.saveAsPDF(new File("Exposure" + ".pdf"), 800, 600);
		
		
	}

}

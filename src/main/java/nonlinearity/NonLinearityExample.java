package nonlinearity;

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

public class NonLinearityExample {
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
		double strike1 =  100;
		double strike2 =  90;
		double maturity1 = deltaT *(numberOfTimeSteps -1);
		
		int quantity1 = 1000;
		int quantity2 = 1000;
		
		//V1 = S0 - K\exp(-rT)
		double exactPrice1 = initialValue - strike1 * Math.exp(-collateralRate * maturity1);
		
		//V2 = K2\exp(-rT) - S0
		double exactPrice2 = strike2 * Math.exp(-collateralRate * maturity1) - initialValue;	
		
		/*
		 * Compute the pathwise exposure
		 */
		RandomVariable[] pathwiseExposure1 = new RandomVariable[numberOfTimeSteps];
		pathwiseExposure1[0] = new RandomVariableFromDoubleArray(timeDiscretization.getTime(0), exactPrice1);
		
		RandomVariable[] pathwiseExposure2 = new RandomVariable[numberOfTimeSteps];
		pathwiseExposure2[0] = new RandomVariableFromDoubleArray(timeDiscretization.getTime(0), exactPrice2);
		
		RandomVariable[] pathwiseExposurePtf = new RandomVariable[numberOfTimeSteps];
				
		double cva1 = 0.0;
		double cva2 = 0.0;
		double cvaPtf = 0.0;
		
		double hazardRateCpty = 0.04;
		double lossGivenDefault = 0.6;
		
		for(int i = 0; i < numberOfTimeSteps; i++) {
			double runningTime = timeDiscretization.getTime(i);
			pathwiseExposure1[i] = (bsModel.getAssetValue(runningTime, 0)
					.mult(Math.exp(-collateralRate * (maturity1 - runningTime)))
					.sub(strike1 * Math.exp(-collateralRate * (maturity1 - runningTime)))).mult(quantity1);
			
			pathwiseExposure2[i] = (bsModel.getAssetValue(runningTime, 0)
					.mult(-Math.exp(-collateralRate * (maturity1 - runningTime)))
					.add(strike2 * Math.exp(-collateralRate * (maturity1 - runningTime)))).mult(quantity2);
			
			pathwiseExposurePtf[i] = pathwiseExposure1[i].add(pathwiseExposure2[i]);
			
		}
		
		for(int i = 1; i < numberOfTimeSteps; i++) {
			double t2 = timeDiscretization.getTime(i);
			double t1 = timeDiscretization.getTime(i-1);
			
			double df2 = Math.exp(-collateralRate * t2);
			double df1 = Math.exp(-collateralRate * t1);
			
			double survProbCpty1 = Math.exp(-hazardRateCpty * t1);
			double survProbCpty2 = Math.exp(-hazardRateCpty * t2);
			double pdCpty = survProbCpty1 - survProbCpty2;
			
			RandomVariable currNe1 = (pathwiseExposure1[i].mult(-1.0)).floor(0.0);
			RandomVariable prevNe1 = (pathwiseExposure1[i-1].mult(-1.0)).floor(0.0);
			double ENE1 =  (currNe1.mult(df2).add(prevNe1.mult(df1))).getAverage();
			
			RandomVariable currNe2 = (pathwiseExposure2[i].mult(-1.0)).floor(0.0);
			RandomVariable prevNe2 = (pathwiseExposure2[i-1].mult(-1.0)).floor(0.0);
			double ENE2 =  (currNe2.mult(df2).add(prevNe2.mult(df1))).getAverage();
			
			RandomVariable currNePtf = (pathwiseExposurePtf[i].mult(-1.0)).floor(0.0);
			RandomVariable prevNePtf = (pathwiseExposurePtf[i-1].mult(-1.0)).floor(0.0);
			double ENEPtf =  (currNePtf.mult(df2).add(prevNePtf.mult(df1))).getAverage();
			
			cva1 += ENE1 * pdCpty;
			cva2 += ENE2 * pdCpty;
			cvaPtf += ENEPtf * pdCpty;
		}
		
		cva1 *= lossGivenDefault * 0.5;
		cva2 *= lossGivenDefault * 0.5;
		cvaPtf *= lossGivenDefault * 0.5;
		
		System.out.println("Stand Alone CVA1: "+ cva1);
		System.out.println("Stand Alone CVA2: " + cva2);
		System.out.println("Ptf CVA: "+ cvaPtf);
		System.out.println("Non-Linearity: "+ (cva1+cva2 - cvaPtf));

		
	    DoubleUnaryOperator PFE1 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesEstimatedExposure;
	    	valuesEstimatedExposure = pathwiseExposure1[index];

	    	double exposureQuantile = valuesEstimatedExposure.getQuantile(0.95);

	    	return exposureQuantile;


	    };
	    
	    DoubleUnaryOperator PFE2 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesEstimatedExposure;
	    	valuesEstimatedExposure = pathwiseExposure2[index];

	    	double exposureQuantile = valuesEstimatedExposure.getQuantile(0.95);

	    	return exposureQuantile;


	    };
	    
	    DoubleUnaryOperator PFEptf = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesEstimatedExposure;
	    	valuesEstimatedExposure = pathwiseExposure1[index].add(pathwiseExposure2[index]);

	    	double exposureQuantile = valuesEstimatedExposure.getQuantile(0.95);

	    	return exposureQuantile;


	    };
	    
	    
	    Plot2D plot = new Plot2D(0.0, 0.99, 900,Arrays.asList(
	        new Named<DoubleUnaryOperator>("PFE1", PFE1),
	        new Named<DoubleUnaryOperator>("PFE2", PFE2),
	        new Named<DoubleUnaryOperator>("PFEptf", PFEptf)));
	    plot.setTitle("95% Potential Future Exposure").setXAxisLabel("Time").setYAxisLabel("Exposure").setIsLegendVisible(true);
	    plot.show();
	    plot.saveAsPDF(new File("PFEImpact" + ".pdf"), 800, 600);
	    
	    
	    DoubleUnaryOperator EPE1 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesPositiveExposure = pathwiseExposure1[index].floor(0.0);

	    	double expectedPositiveExposure   = valuesPositiveExposure.getAverage();

	    	return expectedPositiveExposure;


	    };

	    DoubleUnaryOperator EPE2 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesPositiveExposure = pathwiseExposure2[index].floor(0.0);

	    	double expectedPositiveExposure   = valuesPositiveExposure.getAverage();

	    	return expectedPositiveExposure;


	    };

	    DoubleUnaryOperator EPEptf = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesPositiveExposure = pathwiseExposure1[index].add(pathwiseExposure2[index]).floor(0.0);

	    	double expectedPositiveExposure   = valuesPositiveExposure.getAverage();

	    	return expectedPositiveExposure;


	    };

	    
	    Plot2D plot2 = new Plot2D(0.0, 0.99, 900,Arrays.asList(
	        new Named<DoubleUnaryOperator>("EPE1", EPE1),
	        new Named<DoubleUnaryOperator>("EPE2", EPE2),
	        new Named<DoubleUnaryOperator>("EPEptf", EPEptf)));
	    plot2.setTitle("Expected Positive Exposure").setXAxisLabel("Time").setYAxisLabel("Exposure").setIsLegendVisible(true);
	    plot2.show();
	    plot2.saveAsPDF(new File("EPEImpact" + ".pdf"), 800, 600);
	
	
		    
	    DoubleUnaryOperator ENE1 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesNegativeExposure = pathwiseExposure1[index].cap(0.0);

	    	double expectedNegativeExposure   = valuesNegativeExposure.getAverage();

	    	return expectedNegativeExposure;

	    };
	    
	    DoubleUnaryOperator ENE2 = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesNegativeExposure = pathwiseExposure2[index].cap(0.0);

	    	double expectedNegativeExposure   = valuesNegativeExposure.getAverage();

	    	return expectedNegativeExposure;

	    };
	    
	    DoubleUnaryOperator ENEptf = (time) -> {
	    	int index = timeDiscretization.getTimeIndexNearestLessOrEqual(time);
	    	RandomVariable valuesNegativeExposure = pathwiseExposure1[index].add(pathwiseExposure2[index]).cap(0.0);

	    	double expectedNegativeExposure   = valuesNegativeExposure.getAverage();

	    	return expectedNegativeExposure;

	    };
		    

		    
	    Plot2D plot3 = new Plot2D(0.0, 0.99, 900,Arrays.asList(
	    		new Named<DoubleUnaryOperator>("ENE1", ENE1),
	    		new Named<DoubleUnaryOperator>("ENE2", ENE2),
	    		new Named<DoubleUnaryOperator>("ENEptf", ENEptf)));
	    plot3.setTitle("Expected Negative Exposure").setXAxisLabel("Time").setYAxisLabel("Exposure").setIsLegendVisible(true);
	    plot3.show();
	    plot3.saveAsPDF(new File("ENEImpact" + ".pdf"), 800, 600);


	}

}

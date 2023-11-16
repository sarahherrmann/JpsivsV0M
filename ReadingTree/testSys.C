unsigned int get_bit(int num, int position)
{
	bool bit = num & (1 << position);
	return int(bit);
}

void binary(unsigned int num)
{
	for(int i = 256; i > 0; i = i/2) {
		if(num & i)
			std::cout << "1 ";
		else
			std::cout << "0 ";
	}
	std::cout << std::endl;
}

// enum tailParam {kTailData, kTailMC};
// enum bkgFunction {kPol2Exp, kVWG};
// enum fitRange {kfitRange1, kfitRange2};

enum paramFit {kTailParam, kBkgFunction, kFitRange};

void testSys()
{
  //macro to test bit operation for the systematic uncertainties on the Jpsi fit

  for(unsigned int bitmap = 0; bitmap < 8; bitmap++)
  {
    binary(bitmap);

    unsigned int tailParam = get_bit(bitmap,kTailParam);
    unsigned int bkgFunction = get_bit(bitmap,kBkgFunction);
    unsigned int fitRange = get_bit(bitmap,kFitRange);

    printf("tailParam = %d\n", tailParam);
    printf("bkgFunction = %d\n", bkgFunction);
    printf("fitRange = %d\n", fitRange);

    printf("-----------------\n");


  }

}

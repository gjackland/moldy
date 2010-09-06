import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import java.math.BigDecimal;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

// Args4J imports
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.kohsuke.args4j.spi.BooleanOptionHandler;

public final class Compare
{
	@Option(
		name	= "-d",
		usage	= "Set number of siginficant digits that numbers have to agree to to be considered the same"
	)
	private int		significantDigits	= -1;
	
	@Option(
		name	= "-i",
		usage 	= "Ignore lines that contain the following text"
	)
	private String	ignoreLinesWith		= "";
	
    // Get the list of command line arguments other than options
    @Argument
    private List< String > arguments = new ArrayList< String >();
	
	public static void
	main( final String[] args )
	{
		final Compare comp = new Compare();
		comp.doCompare( args );
	}
	
	public final boolean
	checkFileExists( final File inputFile )
	{
		if( !inputFile.exists() )
		{
			System.out.println( "Error: file " + inputFile.getPath() + " doesn't exist" );
			return false;
		}
		return true;
	}
	
	public void
	doCompare( final String[] args )
	{
        final CmdLineParser parser = new CmdLineParser( this );
        
        // Set the width of the console
        parser.setUsageWidth(80);
        
        try
        {
            // parse the arguments.
            parser.parseArgument( args );

            // after parsing arguments, you should check
            // if enough arguments are given.
            if( arguments.size() != 2 )
                throw new CmdLineException( "Two files for comparison must be supplied" );

        } catch( CmdLineException e )
        {
            // if there's a problem in the command line,
            // you'll get this exception. this will report
            // an error message.
            System.err.println( e.getMessage() );
            System.err.println("java Compare [options...] arguments...");
            // print the list of available options
            parser.printUsage( System.err );
            System.err.println();

            return;
        }
		
		final File file1 = new File( arguments.get( 0 ) );
		final File file2 = new File( arguments.get( 1 ) );
		
		if( !checkFileExists( file1 ) )
		{
			return;
		}
		if( !checkFileExists( file2 ) )
		{
			return;
		}
		
		BufferedReader file1Reader, file2Reader;
		
		
		// Open the files for reading
		try
		{
			file1Reader		= new BufferedReader( new FileReader( file1 ) );
			file2Reader		= new BufferedReader( new FileReader( file2 ) );
		}
		catch( final Exception e )
		{
			System.err.println( "Failed to open file for reading" );
			e.printStackTrace();
			return;
		}
		
		// Read line by line comparing contents
		String	file1Line, file2Line;
		int		lineNo = 0;
		try
		{
			while( file1Reader.ready() && file2Reader.ready() )
			{
				lineNo++;
				
				file1Line = file1Reader.readLine();
				file2Line = file2Reader.readLine();
				
				if( !file1Line.equals( file2Line ) )
				{
					// Check if we're using the ignore flag and either of the files contains the ignore string
					if( !( ignoreLinesWith.length() != 0 &&
						( file1Line.contains( ignoreLinesWith ) || file2Line.contains( ignoreLinesWith ) ) ) )
					{
						final String difference = deepCompare( file1Line, file2Line );
						
						if( difference.length() != 0 )
						{
							System.out.println( lineNo + ":\n" + difference );
						}
					}
				}
			}
		}
		catch( Exception e )
		{
			e.printStackTrace();
			return;
		}
	}
	
	// Compare the tokens in the strings to see how they differ
	private String
	deepCompare( final String s1, final String s2 )
	{
		StringBuffer line1Diff = new StringBuffer(), line2Diff = new StringBuffer();		// To store tokens where the two lines differ
		
		final StringTokenizer st1 = new StringTokenizer( s1 );
		final StringTokenizer st2 = new StringTokenizer( s2 );
		
		String token1, token2;
		while( st1.hasMoreTokens() && st2.hasMoreTokens() )
		{
			token1 = st1.nextToken();
			token2 = st2.nextToken();
			
			if( !token1.equals( token2 ) )
			{
				boolean different = true;
				// If they are numbers check if they differ by more than the tolerance
				try
				{
					final BigDecimal bd1 = new BigDecimal( token1 );
					final BigDecimal bd2 = new BigDecimal( token2 );
					
					// Only check if users has asked to perform significant digits check
					if( significantDigits != -1 )
					{
						different = !withinTolerance( bd1, bd2 );
					}
				}
				catch( NumberFormatException e )
				{
					different = true;
				}
				if( different )
				{
					line1Diff.append( token1 + " " );
					line2Diff.append( token2 + " " );
				}
			}
		}
		
		// Check if one line has more tokens than another
		while( st1.hasMoreTokens() )
		{
			line1Diff.append( st1.nextToken() + " " );
		}
		while( st2.hasMoreTokens() )
		{
			line2Diff.append( st2.nextToken() + " " );
		}
		
		if( line1Diff.length() != 0 || line2Diff.length() != 0 )
		{
			return line1Diff + "\n" + line2Diff;
		}
		
		return "";
	}
	
	private boolean
	withinTolerance( final BigDecimal bd1, final BigDecimal bd2 )
	{
		
		// Figure out the scale of the lower number, the precision minus the scale gives us
		// the exponent had we written the number as 0.xxxE<scaleFactor>
		final int scaleFactor = Math.min( bd1.precision() - bd1.scale(), bd2.precision() - bd2.scale() );
		
		final BigDecimal remainderScaled = bd1.subtract( bd2 ).abs().scaleByPowerOfTen( significantDigits - scaleFactor );
		
		/*
		System.out.println(
			"d1:\t" + bd1 +
			"\nscale:\t" + bd1.scale() +
			"\nprecision:\t" + bd1.precision() +
			"\nd2:\t" + bd2 +
			"\nscale:\t" + bd2.scale() +
			"\nprecision:\t" + bd2.precision() +
			"\nremainder:\t" + remainder +
			"\nremainder scaled:\t" + remainderScaled );
		*/
		
		// Now if the remainder is greater than 1 then we know there is a variance in a significant
		// digit at or before the the one we are checking against
		if( remainderScaled.compareTo( BigDecimal.ONE ) == 1.0 )
		{
			return false;
		}
		return true;
	}
}
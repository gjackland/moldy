
import java.io.*;

import java.net.URI;

import java.lang.NumberFormatException;

import java.util.StringTokenizer;


class FixParams
{
	final static String TMP_PARAMS_APPEND	= ".tmp";
	final static String NUM_ATOMS_TOKEN		= "nm";
	final static String SYSTEM_IN_TOKEN		= "file_system";
	
	public static final boolean
	checkFileExists( File inputFile )
	{
		if( !inputFile.exists() )
		{
			System.out.println( "Error: file " + inputFile.getPath() + " doesn't exist" );
			return false;
		}
		return true;
	}
	
	public static final File
	getRelativePath( final File from, final File to )
	{
		final URI fromURI	= from.getAbsoluteFile().toURI();
		final URI toURI		= to.getAbsoluteFile().toURI();
		final File relativeFile	= new File( fromURI.relativize( toURI ) );
		
		final CharSequence toRemove	= ( CharSequence )( from.getAbsoluteFile().getParent().toString() + File.separator );
		
		final String relativePath = relativeFile.toString().replace( toRemove, ( CharSequence )"" );
		
		return new File( relativePath );
	}
		
	public static void main( String args[] )
	{
		if( args.length != 2 )
		{
			System.out.println( "Usage: java FixParams [system input file] [params input file]"  );
			return;
		}
		
		final File systemIn 	= new File( args[ 0 ] );
		final File paramsIn 	= new File( args[ 1 ] );
		final File paramsInTmp 	= new File( paramsIn.getPath() + TMP_PARAMS_APPEND );
		
		if( !checkFileExists( systemIn ) ) return;
		if( !checkFileExists( paramsIn ) ) return;
		
		
		StringTokenizer	st;
		Integer numAtoms, n_dim;
		
		// First get the require information from the system input file
		try
		{
			final BufferedReader systemInBuffered	= new BufferedReader( new FileReader( systemIn ) );
			st										= new StringTokenizer( systemInBuffered.readLine() );
			final String numAtomsString				= st.nextToken();
			numAtoms								= Integer.valueOf( numAtomsString );
			// Now we need to get the number of times that the cell is repeated in this lattice which should be on the next line
			st										= new StringTokenizer( systemInBuffered.readLine() );
			for( int i = 0; i < 3; ++i )
			{
				n_dim		= Integer.valueOf( st.nextToken() );
				numAtoms	*= n_dim;
			}
		}
		catch( Exception e )
		{
			System.out.println( "Failed to determine number of atoms in " + args[ 0 ] );
			e.printStackTrace();
			return;
		}
		
		BufferedReader paramsInBuffered;
		BufferedWriter paramsInTmpBuffered;
		
		try
		{
			paramsInBuffered	= new BufferedReader( new FileReader( paramsIn ) );
		}
		catch( Exception e )
		{
			System.out.println( "Failed to open " + paramsIn + " for buffered input ");
			e.printStackTrace();
			return;
		}
		
		try
		{
			paramsInTmpBuffered = new BufferedWriter( new FileWriter( paramsInTmp ) );
		}
		catch( Exception e )
		{
			System.out.println( "Failed to open " + paramsIn + TMP_PARAMS_APPEND + " for buffered output");
			e.printStackTrace();
			return;
		}
		
		
		// Now locate the information in the params.in file and make sure the two are in agreement
		String line, firstToken;
		try
		{
			while( paramsInBuffered.ready() )
			{
				try
				{
					line = paramsInBuffered.readLine();
					
					st			= new StringTokenizer( line );
					if( st.hasMoreTokens() )
					{
						firstToken	= st.nextToken();
						
						if( firstToken.equals( NUM_ATOMS_TOKEN ) )
						{
							st.nextToken();		// Skip over the '='
							final String oldNumAtoms = st.nextToken();
							
							// Split the line either side of the old number of atoms
							String oldLine[] = line.split( oldNumAtoms );
							
							line = oldLine[ 0 ] + numAtoms;
							// Add the rest of the line back i.e. a comment
							if( oldLine.length > 1 )
							{
								line += oldLine[ 1 ];
							}
						}
						else if( firstToken.equals( SYSTEM_IN_TOKEN ) )
						{
							st.nextToken();		// Skip over the '='

							final String oldSystemIn = st.nextToken();

							// Have to use the \Q (start) and \E (end) special characters as
							// split matches as a regex and therefore if the token
							// has any metacharacters this will cause havok!						
							// Split the line either side of the old system input file
							String oldLine[] = line.split( "\\Q" + oldSystemIn + "\\E");
							
							final File newSystemIn = getRelativePath( paramsIn, systemIn );
							
							line = oldLine[ 0 ] + "\"" + newSystemIn + "\"";
							// Add the rest of the line back i.e. a comment
							if( oldLine.length > 1 )
							{
								line += oldLine[ 1 ];
							}
							System.out.println( "Changing " + SYSTEM_IN_TOKEN + " from\n" + oldSystemIn + " to\n" + newSystemIn );
						}
					}
					
					// Write the line back out and put a carriage return afterwards
					paramsInTmpBuffered.write( line );		
					paramsInTmpBuffered.newLine();
				}
				catch( Exception e )
				{
					System.out.println( "Failed to read or write line" );
					e.printStackTrace();
					break;
				}
			}
		}
		catch( Exception e )
		{
			System.out.println( "Faied to check if input file is ready to read another line" );
			e.printStackTrace();
		}
	
		try
		{
			paramsInBuffered.close();
			paramsInTmpBuffered.close();
		}
		catch( Exception e)
		{
			System.out.println( "Failed to close file" );
			e.printStackTrace();
			return;
		}
		
		if( paramsIn.delete() )
		{
			if( paramsInTmp.renameTo( paramsIn ) )
			{
				System.out.println( "Sucessfully change number of atoms in " + paramsIn + " to " + numAtoms );
			}
		}
	}
}
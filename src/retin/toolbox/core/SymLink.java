 /**
 * \file SymLink.java
 * \author Philippe H. Gosselin
 * \version 4.0
 */

package retin.toolbox.core;

import java.io.File;

public class SymLink 
{
	public static void create (String file,String dest) throws Exception
	{
/*		File f = new File(file);
		if (f.exists())
			f.delete();*/
		
		String[] command = { "ln" , "-s", "-f", "-n" , dest , file };
		try {
			Runtime.getRuntime().exec(command);
		}
		catch (java.io.IOException e)
		{
			throw new Exception ("Problème lors de l'execution de "+command);
		}
		
	}
	
}

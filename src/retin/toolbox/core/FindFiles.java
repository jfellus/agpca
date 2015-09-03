 /**
 * \file FindFiles.java
 * \author David Gorisse
 * \version 4.0
 */

package retin.toolbox.core;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;


public class FindFiles 
{
//**********************************************************************************************************************
// Look for a file describe by a FilenameFilter in the current directory and all subdirectory
	public static String[] lookInFor(String directory, FilenameFilter filter)
	{
//		System.out.println("main directory : " + directory);
		ArrayList<String> vectFiles = new ArrayList<String>();
		directoryLooping(directory, new File(directory), vectFiles, filter);
		String[] files=new String[vectFiles.size()];
		vectFiles.toArray(files);
		return files;
	}
//**********************************************************************************************************************
// find all files describe by FilenameFilter in the current directory
// files found in the current directory contain the pathname from the main directory
	private static String[] filesOfDir(String mainDir, String currentDir, FilenameFilter filter)
	{

// filter exemple
//			new java.io.FilenameFilter() {
//			public boolean accept(java.io.File dir,String name) {
//				return name.endsWith(".png")
//					|| name.endsWith(".jpg")
//					|| name.endsWith(".jpeg")
//					|| name.endsWith(".gif")
//					|| name.endsWith(".tif")
//					|| name.endsWith(".tiff");
//			}};

		String path=currentDir.substring(mainDir.length());
//	System.out.println("path : " + path);

		if ( !path.endsWith("/") )
			path+="/";
		if ( path.startsWith("/") )
			path=path.substring(1);

//	System.out.println("path : " + path);

		String[] files = new File(currentDir).list(filter);
		
		for (int i=0; i<files.length; i++)
		{
			files[i] = path+files[i];
//			System.out.println("File["+i+"] : " + files[i]);
		}

		return files;
	}
//**********************************************************************************************************************
//Traverse directory Recursively
	private static void directoryLooping(String mainDir, File currentDir, ArrayList<String> vectFiles, FilenameFilter filter)
	{
		if(currentDir.isDirectory())
		{
        	String children[] = currentDir.list();
			for (int i=0; i<children.length; i++)
			{
				File childDir = new File(currentDir,children[i]);
				if (childDir.isDirectory())
					directoryLooping(mainDir, childDir, vectFiles, filter);
			}
	
			String[] files=filesOfDir(mainDir, currentDir.toString(), filter);
			if (files.length != 0 )
				vectFiles.addAll(Arrays.asList(files));
		}
	}

}









	


